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
    double normalizedDist = distance / (size/2.0); // 0 at center, 1 at boundary
    
    if (distance < 1e-10) return 0.0; // Zero force at center
    
    // Force increases as particles approach the boundary
    // Using a quadratic profile: force = k * (r/R)^2 where r is distance from center and R is field radius
    double forceMagnitude = fieldStrength * normalizedDist * normalizedDist;
    
    return forceMagnitude;
}

bool ContainmentField::isParticleContained(const Particle& particle) const {
    std::lock_guard<std::mutex> lock(fieldMutex);
    double x = particle.getX();
    double y = particle.getY();
    double distanceFromCenter = std::sqrt(x*x + y*y);
    double maxRadius = size/2.0;
    
    // Check if particle is within field bounds and has low enough energy to be contained
    return distanceFromCenter < maxRadius && 
           particle.getEnergy() < fieldStrength * (1.0 - distanceFromCenter/maxRadius);
}

void ContainmentField::update(double dt) {
    std::lock_guard<std::mutex> lock(fieldMutex);
    
    // Update field strength with decay
    fieldStrength *= (1.0 - decayRate * dt);
    
    // Update grid data
    double totalEnergy = 0.0;
    for (size_t i = 0; i < fieldData.size(); ++i) {
        fieldData[i] *= (1.0 - decayRate * dt);
        totalEnergy += fieldData[i];
    }
    
    // Update field energy
    fieldEnergy = totalEnergy + fieldStrength * size * size;
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