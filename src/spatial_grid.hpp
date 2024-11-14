#pragma once

#include <vector>
#include <iostream>
#include "cgp/cgp.hpp"
#include "simulation/simulation.hpp"

class spatial_grid_base {
public:
    std::vector<std::vector<particle_element*>> grid;
    float cell_size;
    int grid_width;
    int grid_height;
    cgp::vec3 min;

    spatial_grid_base() = default;

    spatial_grid_base(float cell_size, int width, int height, cgp::vec3 min)
        : cell_size(cell_size), grid_width(width), grid_height(height), min(min) {}

    virtual ~spatial_grid_base() = default;

    virtual int compute_cell_index(const cgp::vec3& position) const = 0;

    void insert_particle(particle_element* particle) {
        int cell_index = compute_cell_index(particle->p);
        if (cell_index < 0 || cell_index >= grid.size()) {
            std::cout << "Invalid cell index: " << cell_index << " for particle " << particle->p << std::endl;
            return;
        }
        particle->cell_index = cell_index;
        grid[cell_index].push_back(particle);
    }

    virtual std::vector<particle_element*> get_neighbors(int cell_index) const = 0;

    void clear_grid() {
        #pragma omp parallel for
        for (int i = 0; i < grid.size(); ++i) {
            grid[i].clear();
        }
    }
};

class spatial_grid_2d : public spatial_grid_base {
public:
    spatial_grid_2d() = default;

    spatial_grid_2d(float cell_size, int width, int height, cgp::vec2 min = {0, 0})
        : spatial_grid_base(cell_size, width, height, cgp::vec3(min.x, min.y, 0)) {
        grid.resize(width * height);
    }

    int compute_cell_index(const cgp::vec3& position) const override {
        int x = (position.x - min.x) / cell_size;
        int y = (position.y - min.y) / cell_size;
        return x + y * grid_width;
    }

    std::vector<particle_element*> get_neighbors(int cell_index) const override {
        std::vector<particle_element*> neighbors;
        int x = cell_index % grid_width;
        int y = cell_index / grid_width;

        for (int dx = -1; dx <= 1; ++dx) {
            for (int dy = -1; dy <= 1; ++dy) {
                int neighbor_x = x + dx;
                int neighbor_y = y + dy;
                if (neighbor_x >= 0 && neighbor_x < grid_width && neighbor_y >= 0 && neighbor_y < grid_height) {
                    int neighbor_index = neighbor_x + neighbor_y * grid_width;
                    neighbors.insert(neighbors.end(), grid[neighbor_index].begin(), grid[neighbor_index].end());
                }
            }
        }
        return neighbors;
    }
};

class spatial_grid_3d : public spatial_grid_base {
public:
    int grid_depth;

    spatial_grid_3d() = default;

    spatial_grid_3d(float cell_size, int width, int height, int depth, cgp::vec3 min = {0, 0, 0})
        : spatial_grid_base(cell_size, width, height, min), grid_depth(depth) {
        grid.resize(width * height * depth);
    }

    int compute_cell_index(const cgp::vec3& position) const override {
        int x = (position.x - min.x) / cell_size;
        int y = (position.y - min.y) / cell_size;
        int z = (position.z - min.z) / cell_size;
        return x + y * grid_width + z * grid_width * grid_height;
    }

    std::vector<particle_element*> get_neighbors(int cell_index) const override {
        std::vector<particle_element*> neighbors;
        int x = cell_index % grid_width;
        int y = (cell_index / grid_width) % grid_height;
        int z = cell_index / (grid_width * grid_height);

        for (int dx = -1; dx <= 1; ++dx) {
            for (int dy = -1; dy <= 1; ++dy) {
                for (int dz = -1; dz <= 1; ++dz) {
                    int neighbor_x = x + dx;
                    int neighbor_y = y + dy;
                    int neighbor_z = z + dz;
                    if (neighbor_x >= 0 && neighbor_x < grid_width &&
                        neighbor_y >= 0 && neighbor_y < grid_height &&
                        neighbor_z >= 0 && neighbor_z < grid_depth) {
                        int neighbor_index = neighbor_x + neighbor_y * grid_width + neighbor_z * grid_width * grid_height;
                        neighbors.insert(neighbors.end(), grid[neighbor_index].begin(), grid[neighbor_index].end());
                    }
                }
            }
        }
        return neighbors;
    }
};