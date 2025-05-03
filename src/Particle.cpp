#include "../include/Particle.h"
#include <cmath>
#include <thread>
#include <chrono>
#include <algorithm>

Particle::Particle(double x, double y, double energy, double radius, double max_energy, double mass)
    : x(x), y(y), vx(0.0), vy(0.0), energy(energy), MAX_ENERGY(max_energy), PARTICLE_RADIUS(radius), mass(mass) {
}

Particle::~Particle() {}

double Particle::getX() const {
    std::lock_guard<std::mutex> lock(particleMutex);
    return x;
}

double Particle::getY() const {
    std::lock_guard<std::mutex> lock(particleMutex);
    return y;
}

void Particle::setPosition(double newX, double newY) {
    std::lock_guard<std::mutex> lock(particleMutex);
    x = newX;
    y = newY;
}

double Particle::getVX() const {
    std::lock_guard<std::mutex> lock(particleMutex);
    return vx;
}

double Particle::getVY() const {
    std::lock_guard<std::mutex> lock(particleMutex);
    return vy;
}

void Particle::setVelocity(double newVX, double newVY) {
    std::lock_guard<std::mutex> lock(particleMutex);
    vx = newVX;
    vy = newVY;
}

double Particle::getEnergy() const {
    std::lock_guard<std::mutex> lock(particleMutex);
    return energy;
}

double Particle::getMaxEnergy() const {
    return MAX_ENERGY;
}

double Particle::getMass() const {
    return mass;
}

void Particle::setEnergy(double newEnergy) {
    std::lock_guard<std::mutex> lock(particleMutex);
    energy = std::clamp(newEnergy, 0.0, MAX_ENERGY);
}

void Particle::addEnergy(double delta) {
    setEnergy(getEnergy() + delta);
}

void Particle::collide(Particle& other) {
    if (this == &other) return;
    std::scoped_lock lock(particleMutex, other.particleMutex);

    // Calculate collision normal
    double dx = other.x - x;
    double dy = other.y - y;
    double dist = std::sqrt(dx*dx + dy*dy);
    if (dist < 1e-10) return;

    double nx = dx / dist;  // Normal vector x
    double ny = dy / dist;  // Normal vector y

    // Relative velocity
    double rvx = other.vx - vx;
    double rvy = other.vy - vy;

    // Relative velocity along normal
    double velAlongNormal = rvx * nx + rvy * ny;
    if (velAlongNormal > 0) return; // Objects moving apart

    // Coefficient of restitution (1.0 for perfectly elastic)
    const double restitution = 1.0;

    // Impulse scalar
    double j = -(1.0 + restitution) * velAlongNormal;
    j /= 1.0/mass + 1.0/other.mass;

    // Apply impulse
    double impulseX = j * nx;
    double impulseY = j * ny;

    vx -= impulseX / mass;
    vy -= impulseY / mass;
    other.vx += impulseX / other.mass;
    other.vy += impulseY / other.mass;

    // Conserve energy by scaling velocities
    double totalEnergyBefore = mass*(vx*vx + vy*vy)/2.0 + other.mass*(other.vx*other.vx + other.vy*other.vy)/2.0;
    double totalEnergyAfter = mass*(vx*vx + vy*vy)/2.0 + other.mass*(other.vx*other.vx + other.vy*other.vy)/2.0;
    
    if (totalEnergyAfter > 0) {
        double scale = std::sqrt(totalEnergyBefore/totalEnergyAfter);
        vx *= scale;
        vy *= scale;
        other.vx *= scale;
        other.vy *= scale;
    }

    // Prevent overlap
    double overlap = PARTICLE_RADIUS + other.PARTICLE_RADIUS - dist;
    if (overlap > 0) {
        double moveX = (overlap/2.0) * nx;
        double moveY = (overlap/2.0) * ny;
        x -= moveX;
        y -= moveY;
        other.x += moveX;
        other.y += moveY;
    }
}

bool Particle::isColliding(const Particle& other) const {
    std::lock_guard<std::mutex> lock1(particleMutex);
    std::lock_guard<std::mutex> lock2(other.particleMutex);
    double dx = x - other.x;
    double dy = y - other.y;
    double distSq = dx*dx + dy*dy;
    double minDist = PARTICLE_RADIUS + other.PARTICLE_RADIUS;
    return distSq < (minDist * minDist);
}
