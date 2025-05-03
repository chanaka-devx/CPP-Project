#pragma once

#include "Particle.h"
#include "ContainmentField.h"
#include "ThreadManager.h"
#include <vector>
#include <memory>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>

struct Config;

class Simulation {
public:
    Simulation(const Config& config);
    ~Simulation();

    void initializeParticles(const Config& config);
    void setContainmentField(std::unique_ptr<ContainmentField> field);

    void start();
    void stop();
    void step();

    void addParticle(std::unique_ptr<Particle> particle);
    void removeEscapedParticles(); 
    size_t getParticleCount() const;
    std::vector<std::unique_ptr<Particle>> getParticlesCopy() const;
    const std::vector<std::unique_ptr<Particle>>& getParticles() const; 

    double getTotalEnergy() const; 

    void setNumThreads(size_t numThreads);
    size_t getNumThreads() const;
    const ThreadManager& getThreadManager() const { return *threadManager; }

    void updatePositions(double timeStep); 
    void applyForces(double timeStep);
    void handleCollisions();

private:
    void workerThread(size_t threadId);  

    struct GridCell {
        std::vector<size_t> particleIndices;
        void clear() { particleIndices.clear(); }
    };
    
    std::vector<GridCell> spatialGrid;
    const double gridCellSize = 2.0;
    size_t gridWidth, gridHeight;
    
    void updateSpatialGrid();
    std::pair<int, int> getGridCoords(double x, double y) const;
    size_t getGridIndex(int gridX, int gridY) const;

    std::vector<std::unique_ptr<Particle>> particles;
    std::unique_ptr<ContainmentField> containmentField;
    std::unique_ptr<ThreadManager> threadManager;
    double fieldSize; 
    const double timeStep; 

    std::vector<std::thread> workerThreads;
    std::mutex simulationMutex;
    std::mutex particleMutex;  
    std::condition_variable cv;
    std::atomic<bool> running{false};
    size_t numThreads;
};