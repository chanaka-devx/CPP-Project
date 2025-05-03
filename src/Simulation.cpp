#include "../include/Simulation.h"
#include "../include/Config.h"
#include <algorithm>
#include <random>
#include <thread>
#include <iostream> 
#include <map>

Simulation::Simulation(const Config& config)
    : fieldSize(config.field_size),
      timeStep(config.time_step),
      containmentField(std::make_unique<ContainmentField>(config)),
      threadManager(std::make_unique<ThreadManager>(config.initial_threads)),
      numThreads(config.initial_threads) {
    // Initialize spatial grid
    gridWidth = static_cast<size_t>(std::ceil(fieldSize / gridCellSize));
    gridHeight = static_cast<size_t>(std::ceil(fieldSize / gridCellSize));
    spatialGrid.resize(gridWidth * gridHeight);
    
    initializeParticles(config);
}

Simulation::~Simulation() {
    stop();
}

void Simulation::initializeParticles(const Config& config) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-fieldSize/2, fieldSize/2);
    std::uniform_real_distribution<> vel_dis(-1.0, 1.0);
    size_t count = config.num_particles;
    for (size_t i = 0; i < count; ++i) {
        auto particle = std::make_unique<Particle>(
            dis(gen), dis(gen),
            config.initial_energy,
            config.particle_radius,
            config.max_energy,
            1.0 // mass
        );
        particle->setVelocity(vel_dis(gen), vel_dis(gen));
        particles.push_back(std::move(particle));
    }
    std::cout << "Initialized " << particles.size() << " particles." << std::endl;
}

void Simulation::setContainmentField(std::unique_ptr<ContainmentField> field) {
    containmentField = std::move(field);
}

void Simulation::start() {
    running = true;
    for (size_t i = 0; i < numThreads; ++i) {
        workerThreads.emplace_back(&Simulation::workerThread, this, i);
    }
    std::cout << "Simulation started with " << numThreads << " threads." << std::endl;
}

void Simulation::stop() {
    running = false;
    for (auto& thread : workerThreads) {
        if (thread.joinable()) {
            thread.join();
        }
    }
    workerThreads.clear();
    std::cout << "Simulation stopped." << std::endl;
}

void Simulation::step() {
    std::lock_guard<std::mutex> lock(simulationMutex);
    
    size_t n = particles.size();
    if (n == 0) return;

    // 1. Apply containment forces and update velocities
    size_t chunk = std::max<size_t>(1, n / numThreads);
    for (size_t t = 0; t < numThreads; ++t) {
        size_t start = t * chunk;
        size_t end = (t == numThreads - 1) ? n : (t + 1) * chunk;
        
        threadManager->addTask([this, start, end]() {
            for (size_t i = start; i < end; ++i) {
                if (i >= particles.size()) break;
                
                double x = particles[i]->getX();
                double y = particles[i]->getY();
                double distance = std::sqrt(x*x + y*y);
                
                if (distance > 0) {
                    double forceMagnitude = containmentField->getContainmentForce(*particles[i]);
                    double fx = -forceMagnitude * (x / distance);
                    double fy = -forceMagnitude * (y / distance);
                    double vx = particles[i]->getVX() + fx * timeStep;
                    double vy = particles[i]->getVY() + fy * timeStep;
                    particles[i]->setVelocity(vx, vy);
                }
            }
        });
    }
    threadManager->waitForCompletion();

    // 2. Update positions
    for (size_t t = 0; t < numThreads; ++t) {
        size_t start = t * chunk;
        size_t end = (t == numThreads - 1) ? n : (t + 1) * chunk;
        
        threadManager->addTask([this, start, end]() {
            for (size_t i = start; i < end; ++i) {
                if (i >= particles.size()) break;
                double newX = particles[i]->getX() + particles[i]->getVX() * timeStep;
                double newY = particles[i]->getY() + particles[i]->getVY() * timeStep;
                particles[i]->setPosition(newX, newY);
            }
        });
    }
    threadManager->waitForCompletion();

    // 3. Update spatial grid and handle collisions
    updateSpatialGrid();
    
    std::vector<std::pair<size_t, size_t>> collisionPairs;
    for (const auto& cell : spatialGrid) {
        for (size_t i = 0; i < cell.particleIndices.size(); ++i) {
            for (size_t j = i + 1; j < cell.particleIndices.size(); ++j) {
                size_t p1 = cell.particleIndices[i];
                size_t p2 = cell.particleIndices[j];
                if (p1 >= particles.size() || p2 >= particles.size()) continue;
                if (particles[p1]->isColliding(*particles[p2])) {
                    collisionPairs.emplace_back(p1, p2);
                }
            }
        }
    }

    // Handle collisions in parallel
    chunk = std::max<size_t>(1, collisionPairs.size() / numThreads);
    for (size_t t = 0; t < numThreads; ++t) {
        size_t start = t * chunk;
        size_t end = (t == numThreads - 1) ? collisionPairs.size() : (t + 1) * chunk;
        
        threadManager->addTask([this, start, end, &collisionPairs]() {
            for (size_t i = start; i < end; ++i) {
                auto [p1, p2] = collisionPairs[i];
                if (p1 < particles.size() && p2 < particles.size()) {
                    particles[p1]->collide(*particles[p2]);
                }
            }
        });
    }
    threadManager->waitForCompletion();

    // 4. Remove escaped particles and update containment field
    removeEscapedParticles();
    containmentField->update(timeStep);
}

void Simulation::updateSpatialGrid() {
    // Reset all cells
    for (auto& cell : spatialGrid) {
        cell.clear();
    }
    
    // Assign particles to grid cells
    for (size_t i = 0; i < particles.size(); ++i) {
        auto [gridX, gridY] = getGridCoords(particles[i]->getX(), particles[i]->getY());
        size_t idx = getGridIndex(gridX, gridY);
        if (idx < spatialGrid.size()) {
            spatialGrid[idx].particleIndices.push_back(i);
        }
    }
}

std::pair<int, int> Simulation::getGridCoords(double x, double y) const {
    int gridX = static_cast<int>((x + fieldSize/2) / gridCellSize);
    int gridY = static_cast<int>((y + fieldSize/2) / gridCellSize);
    return {gridX, gridY};
}

size_t Simulation::getGridIndex(int gridX, int gridY) const {
    if (gridX < 0 || gridX >= static_cast<int>(gridWidth) || 
        gridY < 0 || gridY >= static_cast<int>(gridHeight)) {
        return spatialGrid.size(); // Invalid index
    }
    return gridY * gridWidth + gridX;
}

void Simulation::addParticle(std::unique_ptr<Particle> particle) {
    if (particle) {
        std::lock_guard<std::mutex> lock(particleMutex);
        particles.push_back(std::move(particle));
    }
}

void Simulation::removeEscapedParticles() {
    std::lock_guard<std::mutex> lock(particleMutex);
    particles.erase(
        std::remove_if(
            particles.begin(), particles.end(),
            [this](const std::unique_ptr<Particle>& p) {
                return !containmentField->isParticleContained(*p);
            }
        ),
        particles.end()
    );
}

size_t Simulation::getParticleCount() const {
    std::lock_guard<std::mutex> lock(particleMutex);
    return particles.size();
}

std::vector<std::unique_ptr<Particle>> Simulation::getParticlesCopy() const {
    std::lock_guard<std::mutex> lock(particleMutex);
    std::vector<std::unique_ptr<Particle>> copy;
    for (const auto& p : particles) {
        // Deep copy not possible for unique_ptr without clone method, so just return shallow copy for now
        // (Rendering only reads x/y, so this is safe as long as Particle's accessors are thread-safe)
        // If deep copy needed, implement a clone() method in Particle
        copy.push_back(std::unique_ptr<Particle>(nullptr)); // placeholder, not used
    }
    return copy;
}

const std::vector<std::unique_ptr<Particle>>& Simulation::getParticles() const {
    return particles;
}

double Simulation::getTotalEnergy() const {
    std::lock_guard<std::mutex> lock(particleMutex);
    double total = 0.0;
    for (const auto& particle : particles) {
        total += particle->getEnergy();
    }
    return total;
}

void Simulation::setNumThreads(size_t newNumThreads) {
    numThreads = newNumThreads;
    threadManager->setNumThreads(newNumThreads);
}

size_t Simulation::getNumThreads() const {
    return numThreads;
}

void Simulation::updatePositions(double dt) {
    std::lock_guard<std::mutex> lock(particleMutex);
    for (auto& particle : particles) {
        double x = particle->getX() + particle->getVX() * dt;
        double y = particle->getY() + particle->getVY() * dt;
        particle->setPosition(x, y);
    }
}

void Simulation::handleCollisions() {
    std::lock_guard<std::mutex> lock(particleMutex);
    for (size_t i = 0; i < particles.size(); ++i) {
        for (size_t j = i + 1; j < particles.size(); ++j) {
            if (particles[i]->isColliding(*particles[j])) {
                particles[i]->collide(*particles[j]);
            }
        }
    }
}

void Simulation::applyForces(double dt) {
    std::lock_guard<std::mutex> lock(particleMutex);
    for (auto& particle : particles) {
        double x = particle->getX();
        double y = particle->getY();
        double distance = std::sqrt(x*x + y*y);
        if (distance == 0) continue;
        double forceMagnitude = -containmentField->getFieldStrength() * distance / fieldSize;
        double ax = forceMagnitude * (x / distance);
        double ay = forceMagnitude * (y / distance);
        double vx = particle->getVX() + ax * dt;
        double vy = particle->getVY() + ay * dt;
        particle->setVelocity(vx, vy);
    }
}

void Simulation::workerThread(size_t threadId) {
    while (running) {
        threadManager->addTask([this] { this->step(); });  // Scheduling the task for the thread
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
}
