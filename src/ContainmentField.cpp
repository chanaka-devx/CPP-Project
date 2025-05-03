#include "../include/ContainmentField.h"
#include "../include/Particle.h"
#include "../include/Config.h"
#include <cmath>
#include <algorithm>

ContainmentField::ContainmentField(const Config& config)
    : size(config.field_size), fieldStrength(config.initial_strength), decayRate(config.initial_decay_rate), GRID_SIZE(config.field_grid_size), fieldEnergy(0.0) {
    initializeField();
}

ContainmentField::~ContainmentField() {
    
}

void ContainmentField::initializeField() {
    fieldData.resize(GRID_SIZE * GRID_SIZE, 0.0); 
}

double ContainmentField::getContainmentForce(const Particle& particle) const {
    double x = particle.getX();
    double y = particle.getY();
    
    // Calculate distance from center as a fraction of field size
    double distance = std::sqrt(x*x + y*y);
    double maxRadius = size/2.0;
    double normalizedDist = distance / maxRadius; // 0 at center, 1 at boundary
    
    if (distance < 1e-10) return 0.0; // Zero force at center
    
    // Force increases exponentially as particles approach the boundary
    // Using profile: force = k * exp(alpha * (r/R - 1))
    // This gives a steep increase near the boundary while maintaining smooth behavior
    const double alpha = 2.0; // Controls how sharply the force increases
    double forceMagnitude = fieldStrength * std::exp(alpha * (normalizedDist - 1.0));
    
    // Scale force by particle energy - higher energy particles experience stronger containment
    forceMagnitude *= (1.0 + particle.getEnergy() / particle.getMaxEnergy());
    
    return forceMagnitude;
}

bool ContainmentField::isParticleContained(const Particle& particle) const {
    std::lock_guard<std::mutex> lock(fieldMutex);
    double x = particle.getX();
    double y = particle.getY();
    double distance = std::sqrt(x*x + y*y);
    double maxRadius = size/2.0;
    
    // Check if particle is within field bounds
    if (distance >= maxRadius) return false;
    
    // Check if particle's energy is low enough to be contained
    double normalizedDist = distance / maxRadius;
    double containmentThreshold = fieldStrength * (1.0 - normalizedDist * normalizedDist);
    return particle.getEnergy() < containmentThreshold;
}

void ContainmentField::update(double dt) {
    std::lock_guard<std::mutex> lock(fieldMutex);
    
    // Update field strength with decay
    double oldStrength = fieldStrength;
    fieldStrength *= (1.0 - decayRate * dt);
    
    // Calculate energy lost from field decay
    double energyLost = (oldStrength - fieldStrength) * size * size;
    
    // Update grid data and accumulate total field energy
    double totalGridEnergy = 0.0;
    for (size_t i = 0; i < fieldData.size(); ++i) {
        fieldData[i] *= (1.0 - decayRate * dt);
        totalGridEnergy += fieldData[i];
    }
    
    // Total field energy is sum of uniform field energy and local perturbations
    fieldEnergy = (fieldStrength * size * size) + totalGridEnergy;
}

void ContainmentField::setFieldStrength(double strength) {
    std::lock_guard<std::mutex> lock(fieldMutex);
    fieldStrength = strength;
}

double ContainmentField::getFieldStrength() const {
    std::lock_guard<std::mutex> lock(fieldMutex);
    return fieldStrength;
}

void ContainmentField::setDecayRate(double rate) {
    std::lock_guard<std::mutex> lock(fieldMutex);
    decayRate = rate;
}

double ContainmentField::getDecayRate() const {
    std::lock_guard<std::mutex> lock(fieldMutex);
    return decayRate;
}

double ContainmentField::getSize() const {
    return size;
}

double ContainmentField::getFieldEnergy() const {
    std::lock_guard<std::mutex> lock(fieldMutex);
    return fieldEnergy;
}