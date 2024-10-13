#include "grid_container.hpp"

#include <functional>

int coords_to_id(cgp::vec3 coords, cgp::vec3 size)
{
    return coords.x + (size.x * coords.y) + (size.x * size.y * coords.z);
}

block_container::block_container(cgp::vec3 id, cgp::vec3 size) :
    particles_(),
    id_(id),
    neighbours_(27)
{
    int i = 0;
    for (int x = -1; x <= 1; x++) {
        for (int y = -1; y <= 1; y++) {
            for (int z = -1; z <= 1; z++) {
                cgp::vec3 neighbour = id + cgp::vec3(x, y, z);
                if (neighbour.x < 0 || neighbour.x >= size.x || neighbour.y < 0 || neighbour.y >= size.y || neighbour.z < 0 || neighbour.z >= size.z) {
                    neighbours_.at(i) = -1;
                }
                else {
                    neighbours_.at(i) = coords_to_id(id + cgp::vec3(x, y, z), size);
                }
                i++;
            }
        }
    }
}

cgp::numarray<int> block_container::get_neighbours_blocks()
{
    return neighbours_;
}

void block_container::add_particle(std::shared_ptr<particle_element> particle)
{
    particles_.push_back(particle);
}

grid_container::grid_container(cgp::vec3 volume_size, cgp::vec3 p, float h)
    : p_(p)
    , h_(h)
{
    size = volume_size / (h * 2.0);
    size = cgp::vec3(std::floor(size.x), std::floor(size.y), std::floor(size.z));
    size = size + cgp::vec3(1, 1, 1);

    blocks.clear();
    for (int x = 0; x < size.x; x++) {
        for (int y = 0; y < size.y; y++) {
            for (int z = 0; z < size.z; z++) {
                blocks.push_back(std::make_shared<block_container>(cgp::vec3(x, y, z), size));
            }
        }
    }
}

void grid_container::clear()
{
    for (int i = 0; i < blocks.size(); i++) {
        blocks.at(i)->clear();
    }
}

cgp::numarray<int> grid_container::get_neighbours_blocks(std::shared_ptr<block_container> block)
{
    int id = coords_to_id(block->get_id(), size);
    return blocks.at(id)->get_neighbours_blocks();
}

void grid_container::add_particle(std::shared_ptr<particle_element> particle)
{
    cgp::vec3 coords = particle->p;
    cgp::vec3 id = (coords - p_) / (h_ * 2.0);
    id = cgp::vec3(std::floor(id.x), std::floor(id.y), std::floor(id.z));
    int block_id = coords_to_id(id, size);

    particle->block = blocks.at(block_id);
    blocks.at(block_id)->add_particle(particle);
}