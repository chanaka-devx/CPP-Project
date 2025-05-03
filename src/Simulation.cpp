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

void Simulation::calculateForces() {
    // Resize force vector if needed
    if (particleForces.size() != particles.size()) {
        particleForces.resize(particles.size());
    }
    
    // Reset forces
    for (auto& force : particleForces) {
        force.fx = 0.0;
        force.fy = 0.0;
    }
    
    // Calculate forces in parallel
    size_t n = particles.size();
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
                    // Force points toward center (negative position components)
                    double fx = -forceMagnitude * (x / distance);
                    double fy = -forceMagnitude * (y / distance);
                    
                    // Atomically accumulate forces
                    particleForces[i].fx.store(fx, std::memory_order_relaxed);
                    particleForces[i].fy.store(fy, std::memory_order_relaxed);
                }
            }
        });
    }
    threadManager->waitForCompletion();
}

void Simulation::applyForces() {
    size_t n = particles.size();
    size_t chunk = std::max<size_t>(1, n / numThreads);
    
    for (size_t t = 0; t < numThreads; ++t) {
        size_t start = t * chunk;
        size_t end = (t == numThreads - 1) ? n : (t + 1) * chunk;
        
        threadManager->addTask([this, start, end]() {
            for (size_t i = start; i < end; ++i) {
                if (i >= particles.size()) break;
                
                double fx = particleForces[i].fx.load(std::memory_order_relaxed);
                double fy = particleForces[i].fy.load(std::memory_order_relaxed);
                
                // Update velocities based on forces
                double vx = particles[i]->getVX() + fx * timeStep;
                double vy = particles[i]->getVY() + fy * timeStep;
                particles[i]->setVelocity(vx, vy);
            }
        });
    }
    threadManager->waitForCompletion();
}

void Simulation::step() {
    std::lock_guard<std::mutex> lock(simulationMutex);
    
    size_t n = particles.size();
    if (n == 0) return;

    // 1. Calculate and apply forces
    calculateForces();
    applyForces();

    // 2. Update positions
    size_t chunk = std::max<size_t>(1, n / numThreads);
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
    // Check each grid cell and its neighbors
    for (size_t y = 0; y < gridHeight; ++y) {
        for (size_t x = 0; x < gridWidth; ++x) {
            size_t cellIdx = y * gridWidth + x;
            const auto& cell = spatialGrid[cellIdx];
            
            // Check within the same cell
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
            
            // Check neighboring cells (right and bottom neighbors only to avoid duplicates)
            const int dx[] = {1, 0, 1};  // right, bottom, bottom-right
            const int dy[] = {0, 1, 1};
            
            for (int dir = 0; dir < 3; ++dir) {
                int nx = static_cast<int>(x) + dx[dir];
                int ny = static_cast<int>(y) + dy[dir];
                
                if (nx >= 0 && nx < static_cast<int>(gridWidth) && 
                    ny >= 0 && ny < static_cast<int>(gridHeight)) {
                    size_t neighborIdx = ny * gridWidth + nx;
                    const auto& neighborCell = spatialGrid[neighborIdx];
                    
                    // Check collisions between current cell and neighbor cell
                    for (size_t p1 : cell.particleIndices) {
                        for (size_t p2 : neighborCell.particleIndices) {
                            if (p1 >= particles.size() || p2 >= particles.size()) continue;
                            if (particles[p1]->isColliding(*particles[p2])) {
                                collisionPairs.emplace_back(p1, p2);
                            }
                        }
                    }
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

void Simulation::recalculateGridSize() {
    // Use avg particle radius * 2 as the cell size for efficient collision detection
    if (particles.empty()) {
        gridCellSize = 2.0; // Default size
        return;
    }

    double avgRadius = 0.0;
    for (const auto& p : particles) {
        avgRadius += p->PARTICLE_RADIUS;
    }
    avgRadius /= particles.size();
    
    // Set grid cell size to be slightly larger than two particle diameters
    gridCellSize = avgRadius * 4.0;
    
    // Update grid dimensions
    gridWidth = static_cast<size_t>(std::ceil(fieldSize / gridCellSize));
    gridHeight = static_cast<size_t>(std::ceil(fieldSize / gridCellSize));
    
    // Resize grid
    spatialGrid.resize(gridWidth * gridHeight);
}

void Simulation::updateSpatialGrid() {
    // Recalculate grid size periodically or when particle count changes significantly
    static size_t lastParticleCount = 0;
    if (std::abs(static_cast<long>(particles.size()) - static_cast<long>(lastParticleCount)) > particles.size() / 10) {
        recalculateGridSize();
        lastParticleCount = particles.size();
    }
    
    // Clear all cells
    for (auto& cell : spatialGrid) {
        cell.clear();
    }
    
    // Reserve space in cells based on average density
    size_t avgParticlesPerCell = (particles.size() / spatialGrid.size()) + 1;
    for (auto& cell : spatialGrid) {
        cell.particleIndices.reserve(avgParticlesPerCell * 2);
    }
    
    // Distribute particles to grid cells in parallel
    size_t n = particles.size();
    size_t chunk = std::max<size_t>(1, n / numThreads);
    
    std::vector<std::vector<std::pair<size_t, size_t>>> threadLocalAssignments(numThreads);
    for (size_t t = 0; t < numThreads; ++t) {
        size_t start = t * chunk;
        size_t end = (t == numThreads - 1) ? n : (t + 1) * chunk;
        
        threadManager->addTask([this, start, end, t, &threadLocalAssignments]() {
            auto& assignments = threadLocalAssignments[t];
            assignments.reserve(end - start);
            
            for (size_t i = start; i < end; ++i) {
                if (i >= particles.size()) break;
                
                auto [gridX, gridY] = getGridCoords(particles[i]->getX(), particles[i]->getY());
                size_t idx = getGridIndex(gridX, gridY);
                if (idx < spatialGrid.size()) {
                    assignments.emplace_back(idx, i);
                }
            }
        });
    }
    threadManager->waitForCompletion();
    
    // Merge thread-local assignments into the spatial grid
    for (const auto& assignments : threadLocalAssignments) {
        for (const auto& [cellIdx, particleIdx] : assignments) {
            spatialGrid[cellIdx].particleIndices.push_back(particleIdx);
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
    copy.reserve(particles.size());
    for (const auto& p : particles) {
        // Create a lightweight copy with just position data for rendering
        copy.push_back(std::make_unique<Particle>(
            p->getX(), p->getY(),
            p->getEnergy(),
            p->PARTICLE_RADIUS,
            p->getMaxEnergy()
        ));
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
