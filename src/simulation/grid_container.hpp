#pragma once

#include "cgp/cgp.hpp"

#include "particle_element.hpp"

class block_container
{
private:
    cgp::numarray<std::shared_ptr<particle_element>> particles_;
    cgp::vec3 id_;

    cgp::numarray<int> neighbours_;

public:
    block_container() = default;
    block_container(cgp::vec3 id, cgp::vec3 size);

    int size() { return particles_.size(); }
    std::shared_ptr<particle_element> at(int id) { return particles_.at(id); }
    cgp::vec3 get_id() { return id_; }

    void clear() { particles_.clear(); }
    cgp::numarray<int> get_neighbours_blocks();

    void add_particle(std::shared_ptr<particle_element> particle);
};

class grid_container
{
private:
    cgp::numarray<std::shared_ptr<block_container>> blocks;
    cgp::vec3 size;
    cgp::vec3 p_;
    float h_;

public:
    grid_container() = default;
    grid_container(cgp::vec3 volume_size, cgp::vec3 p, float h);

    void clear();

    cgp::numarray<int> get_neighbours_blocks(std::shared_ptr<block_container> block);

    std::shared_ptr<block_container> at(int id) { return blocks.at(id); }

    void add_particle(std::shared_ptr<particle_element> particle);
};