#pragma once

#include <map>
#include <list>

#include "cgp/cgp.hpp"
#include "environment.hpp"
#include "utils.hpp"
#include "camera_controller.hpp"
#include "spatial_grid.hpp"

#include "simulation/simulation.hpp"
#include "simulation/simulation_3d.hpp"
#include "simulation/simulation_2d.hpp"

using cgp::mesh_drawable;

struct gui_parameters {
	int right_click_action = SPAWN_PARTICLES;

	int spawn_particle_number = 10;
	float spawn_particle_radius = 0.3f;
	int spawn_particle_type = 0;

	bool new_fluid_creation = false;
	char new_fluid_name[64] = "New fluid";
	float new_fluid_mass = 1.0f;
	float new_fluid_viscosity = 0.01f;
	vec3 new_fluid_color = { 1.0, 1.0, 1.0 };
	//std::vector<std::shared_ptr<bool>> new_fluid_soluble;

	float force_strength = 1.0f;

	bool display_color = true;
	bool display_particles = true;
	bool display_radius = false;

	int grid_size = 50;

	int color_type = FLUID_COLOR;
	float color_smoothing_radius = 0.12f;

	float threshold_min = 0.0f;
	float threshold_max = 1.0f;
	vec3 color_min = { 0,0,1 };
	vec3 color_max = { 1,0,0 };
};

// The structure of the custom scene
struct scene_structure : cgp::scene_inputs_generic {

	int dimension = DIM_2D;
	
	// ****************************** //
	// Elements and shapes of the scene
	// ****************************** //
	camera_controller camera_control;
	camera_projection cam_projection;
	window_structure window;

	mesh_drawable global_frame;          // The standard global frame
	environment_structure environment;   // Standard environment controler
	input_devices inputs;                // Storage for inputs status (mouse, keyboard, window dimension)
	gui_parameters gui;                  // Standard GUI element storage
	
	// ****************************** //
	// Elements and shapes of the scene
	// ****************************** //
	cgp::timer_basic timer;

	sph_parameters_structure sph_parameters; // Physical parameter related to SPH
	cgp::numarray<particle_element> particles;      // Storage of the particles
	spatial_grid grid; // Spatial grid used to accelerate the computation of the SPH
	cgp::numarray<std::shared_ptr<fluid_class>> fluid_classes;       // Storage of the different fluids present in the scene
	cgp::mesh_drawable sphere_particle; // Sphere used to display a particle
	cgp::curve_drawable curve_visual;   // Circle used to display the radius h of influence

	cgp::grid_2D<cgp::vec3> field;      // grid used to represent the volume of the fluid under the particles
	cgp::mesh_drawable field_quad; // quad used to display this field color


	// ****************************** //
	// Functions
	// ****************************** //

	void initialize();    // Standard initialization to be called before the animation loop
	void display_frame(); // The frame display to be called within the animation loop
	void display_gui();   // The display of the GUI, also called within the animation loop

	void new_fluid_gui();

	void update_field_closest(int Nf);
	void update_field_mean(int Nf);
	void update_field_color();
	vec3 get_particle_color(particle_element const& particle);

	void spawn_particle(vec3 const& pos, int fluid_type);
	void spawn_random_type_particle(vec3 const &center);
	void spawn_particles_in_disk(vec3 const& center, float radius, int N, int fluid_type);
	void spawn_particles_in_sphere(vec3 const& center, float radius, int N, int fluid_type);
	void initialize_fluid_classes();
	void initialize_sph();

	void delete_particles_in_disk(vec3 const &center, float radius);
	void add_radial_force(vec3 const &center, float radius);

	void mouse_move_event();
	void mouse_click_event();
	void keyboard_event();
	void idle_frame();

};