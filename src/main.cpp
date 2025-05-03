#include "../include/Simulation.h"
#include "../include/Config.h"
#include <iostream>
#include <chrono>
#include <thread>
#include <fstream>
#include <stdexcept>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

Config loadConfig(const std::string& filename) {
    std::ifstream configFile(filename);
    if (!configFile.is_open()) {
        throw std::runtime_error("Could not open configuration file: " + filename);
    }

    json j;
    try {
        configFile >> j;
    } catch (const json::parse_error& e) {
        throw std::runtime_error("Failed to parse configuration file: " + std::string(e.what()));
    }

    return Config::fromJson(j);
}

void renderASCII(const std::vector<std::unique_ptr<Particle>>& particles, 
                double fieldSize, const Config& cfg) {
    std::vector<std::vector<int>> gridCounts(cfg.grid_height, 
                                           std::vector<int>(cfg.grid_width, 0));

    // Count particles in each grid cell
    for (const auto& particle : particles) {
        if (!particle) continue; // Skip null particles in copy
        double x = particle->getX();
        double y = particle->getY();

        // Convert particle position to grid coordinates
        int col = static_cast<int>((x + fieldSize/2) * cfg.grid_width / fieldSize);
        int row = static_cast<int>((y + fieldSize/2) * cfg.grid_height / fieldSize);

        // Clamp to grid boundaries
        col = std::clamp(col, 0, cfg.grid_width - 1);
        row = std::clamp(row, 0, cfg.grid_height - 1);

        gridCounts[row][col]++;
    }

    // Clear screen and move cursor to top-left
    std::cout << "\033[2J\033[H"; 

    // Draw top border
    std::cout << '+' << std::string(cfg.grid_width, '-') << "+\n";

    // Draw grid
    for (int i = 0; i < cfg.grid_height; ++i) {
        std::cout << '|';
        for (int j = 0; j < cfg.grid_width; ++j) {
            int count = gridCounts[i][j];
            if (count == 0) {
                std::cout << ' ';
            } else {
                int level = std::min(count, cfg.max_density_level);
                auto it = cfg.density_map.find(level);
                std::cout << (it != cfg.density_map.end() ? it->second : ' ');
            }
        }
        std::cout << "|\n";
    }

    // Draw bottom border
    std::cout << '+' << std::string(cfg.grid_width, '-') << "+\n";
    std::cout << std::flush;
}

int main() {
    try {
        const std::string configFilename = "config.json";
        Config config = loadConfig(configFilename);
        std::cout << "Configuration loaded from " << configFilename << std::endl;

        Simulation simulation(config);
        simulation.start();

        const double FRAME_TIME = 1.0 / config.target_fps;
        size_t frameCount = 0;
        auto lastStatTime = std::chrono::high_resolution_clock::now();

        while (simulation.getParticleCount() > 0) {
            auto frameStart = std::chrono::high_resolution_clock::now();

            simulation.step();
            auto particles = simulation.getParticlesCopy();
            renderASCII(particles, config.field_size, config);

            // Calculate and maintain target frame rate
            auto frameEnd = std::chrono::high_resolution_clock::now();
            auto frameDuration = std::chrono::duration<double>(frameEnd - frameStart).count();
            if (frameDuration < FRAME_TIME) {
                std::this_thread::sleep_for(
                    std::chrono::duration<double>(FRAME_TIME - frameDuration)
                );
            }

            // Display statistics every 30 frames
            if (++frameCount % 30 == 0) {
                auto now = std::chrono::high_resolution_clock::now();
                auto elapsed = std::chrono::duration<double>(now - lastStatTime).count();
                double actualFps = (elapsed > 1e-6) ? (30.0 / elapsed) : 0.0;
                lastStatTime = now;

                std::cout << "\nParticles: " << simulation.getParticleCount()
                          << " | Energy: " << simulation.getTotalEnergy()
                          << " | FPS: " << actualFps 
                          << " | Field Strength: " << simulation.getThreadManager().getNumThreads()
                          << std::endl;
            }
        }

        simulation.stop();
        std::cout << "Simulation ended. All particles escaped.\n";

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}