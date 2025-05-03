#include "../include/Particle.h"
#include <cmath>
#include <thread>
#include <chrono>
#include <algorithm>

Particle::Particle(double x, double y, double energy, double radius, double max_energy)
    : x(x), y(y), vx(0.0), vy(0.0), energy(energy), MAX_ENERGY(max_energy), PARTICLE_RADIUS(radius) {
    // Use the provided energy value
}

Particle::~Particle() {
}

double Particle::getX() const {
    return x;
}

double Particle::getY() const {
    return y;
}

void Particle::setPosition(double newX, double newY) {
    x = newX;
    y = newY;
}

double Particle::getVX() const {
    return vx;
}

double Particle::getVY() const {
    return vy;
}

void Particle::setVelocity(double newVX, double newVY) {
    std::lock_guard<std::mutex> lock(particleMutex);
    vx = newVX;
    vy = newVY;
}

double Particle::getEnergy() const {
    return energy;
}

double Particle::getMaxEnergy() const {
    return MAX_ENERGY;
}

void Particle::setEnergy(double newEnergy) {
    energy = std::clamp(newEnergy, 0.0, MAX_ENERGY);
}

void Particle::addEnergy(double delta) {
    setEnergy(energy + delta);
}

void Particle::collide(Particle& other) {
    // Simple elastic collision: swap velocities
    std::lock_guard<std::mutex> lock1(particleMutex);
    std::lock_guard<std::mutex> lock2(other.particleMutex);
    std::swap(vx, other.vx);
    std::swap(vy, other.vy);
}

bool Particle::isColliding(const Particle& other) const {
    double dx = x - other.x;
    double dy = y - other.y;
    double distSq = dx*dx + dy*dy;
    double minDist = PARTICLE_RADIUS + other.PARTICLE_RADIUS;
    return distSq < (minDist * minDist);
}
