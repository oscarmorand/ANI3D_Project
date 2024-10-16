#pragma once

#include <vector>

#include "cgp/cgp.hpp"
#include "simulation/simulation.hpp"

class spatial_grid {
public:
    std::vector<std::vector<particle_element*>> grid;
    float cell_size;
    int grid_width;  // Largeur de la grille en nombre de cellules
    int grid_height;  // Hauteur de la grille en nombre de cellules

    cgp::vec2 min;

    spatial_grid() = default;

    spatial_grid(float cell_size, int width, int height, cgp::vec2 min = {0, 0}) 
        : cell_size(cell_size)
        , grid_width(width)
        , grid_height(height)
        , min(min)
    {
        grid.resize(width * height);  // Initialiser la grille avec la taille totale
    }

    int compute_cell_index(const cgp::vec3& position) {
        int x = (position.x - min.x) / cell_size;
        int y = (position.y - min.y) / cell_size;
        return x + y * grid_width;
    }

    void insert_particle(particle_element* particle) {
        int cell_index = compute_cell_index(particle->p);
        if (cell_index < 0 || cell_index >= grid.size()) {
            std::cout << "Invalid cell index: " << cell_index << " for particle " << particle->p << std::endl;
            return;
        }
        particle->cell_index = cell_index;
        grid[cell_index].push_back(particle);
    }

    std::vector<particle_element*> get_neighbors(int cell_index) {
        std::vector<particle_element*> neighbors;
        int x = cell_index % grid_width;
        int y = cell_index / grid_height;
        // Rechercher dans les cellules voisines (3x3)
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

    void clear_grid() {
        // Pour chaque cellule de la grille, on vide son contenu
        #pragma omp parallel for
        for (int i = 0; i < grid.size(); ++i) {
            grid[i].clear();  // Vide les particules de chaque cellule
        }
    }
};