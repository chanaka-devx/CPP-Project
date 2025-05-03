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
    // Elastic collision: conserve momentum and energy
    double v1x = vx, v1y = vy, v2x = other.vx, v2y = other.vy;
    double m1 = mass, m2 = other.mass;
    vx = (v1x * (m1 - m2) + 2 * m2 * v2x) / (m1 + m2);
    vy = (v1y * (m1 - m2) + 2 * m2 * v2y) / (m1 + m2);
    other.vx = (v2x * (m2 - m1) + 2 * m1 * v1x) / (m1 + m2);
    other.vy = (v2y * (m2 - m1) + 2 * m1 * v1y) / (m1 + m2);
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
