#include "../include/Config.h"
#include <nlohmann/json.hpp>

using json = nlohmann::json;

Config Config::fromJson(const json& j) {
    Config cfg;

    auto get_nested = [&](const json& root, const std::string& group, const std::string& key) {
        if (!root.contains(group) || !root.at(group).contains(key)) {
            throw std::runtime_error("Missing required configuration key: " + group + "." + key);
        }
        return root.at(group).at(key);
    };

    cfg.num_particles = get_nested(j, "simulation", "num_particles");
    cfg.field_size = get_nested(j, "simulation", "field_size");
    cfg.initial_threads = get_nested(j, "simulation", "initial_threads");
    cfg.time_step = get_nested(j, "simulation", "time_step");

    cfg.initial_energy = get_nested(j, "particle", "initial_energy");
    cfg.max_energy = get_nested(j, "particle", "max_energy");
    cfg.particle_radius = get_nested(j, "particle", "radius");

    cfg.initial_strength = get_nested(j, "containment_field", "initial_strength");
    cfg.initial_decay_rate = get_nested(j, "containment_field", "initial_decay_rate");
    cfg.field_grid_size = get_nested(j, "containment_field", "grid_size");

    cfg.target_fps = get_nested(j, "rendering", "target_fps");
    cfg.grid_width = get_nested(j, "rendering", "grid_width");
    cfg.grid_height = get_nested(j, "rendering", "grid_height");
    cfg.max_density_level = get_nested(j, "rendering", "max_density_level");

    const auto& density_map_json = get_nested(j, "rendering", "density_map");
    if (!density_map_json.is_object()) {
        throw std::runtime_error("Configuration error: rendering.density_map must be an object.");
    }

    cfg.density_map.clear();
    for (const auto& [key_str, val] : density_map_json.items()) {
        try {
            int key = std::stoi(key_str);
            if (!val.is_string() || val.get<std::string>().length() != 1) {
                throw std::runtime_error("Configuration error: density_map values must be single characters.");
            }
            cfg.density_map[key] = val.get<std::string>()[0];
        } catch (const std::invalid_argument&) {
            throw std::runtime_error("Configuration error: density_map keys must be integers: " + key_str);
        }
    }

    return cfg;
}