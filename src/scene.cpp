#include "scene.hpp"

using namespace cgp;

void scene_structure::initialize()
{
	if (dimension == DIM_2D)
		cam_projection = cgp::camera_projection::build_orthographic(-1.1f, 1.1f, -1.1f, 1.1f, -10, 10);
	else
		cam_projection = cgp::camera_projection::build_perspective(50.0f * 3.14f / 180, 1.0f, 0.1f, 100.0f);
	cam_projection.update_aspect_ratio(window.aspect_ratio());

	camera_control.initialize(inputs, window); // Give access to the inputs and window global state to the camera controler
	global_frame.initialize_data_on_gpu(mesh_primitive_frame());

	field.resize(50, 50);
	field_quad.initialize_data_on_gpu(mesh_primitive_quadrangle({-1, -1, 0}, {1, -1, 0}, {1, 1, 0}, {-1, 1, 0}));
	field_quad.material.phong = {1, 0, 0};
	field_quad.texture.initialize_texture_2d_on_gpu(field);

	initialize_fluid_classes();
	initialize_sph();
	sphere_particle.initialize_data_on_gpu(mesh_primitive_sphere(1.0, {0, 0, 0}, 10, 10));
	sphere_particle.model.scaling = 0.01f;
	curve_visual.color = {1, 0, 0};
	curve_visual.initialize_data_on_gpu(curve_primitive_circle());
}

void scene_structure::spawn_particle(vec3 const &pos, int fluid_type)
{
	particle_element particle;
	particle.p = pos;
	particle.v = {0, 0, 0};
	particle.f = {0, 0, 0};

	auto const& fluid_class = fluid_classes[fluid_type];
	particle.m = fluid_class->base_m;
	particle.nu = fluid_class->base_nu;
	particle.color = fluid_class->base_color;

	particle.fluid_type = fluid_class;

	particles.push_back(particle);
}

void scene_structure::spawn_random_type_particle(vec3 const &center) {
	float m = rand_uniform();
	if (m == 1.0f)
		m = 0.0f;

	int n = fluid_classes.size();
	if (n == 0)
		return;

	int k = int(m * n);

	spawn_particle(center, k);
}

void scene_structure::spawn_particles_in_disk(vec3 const &center, float radius, int N, int fluid_type)
{
	for (int k = 0; k < N; ++k)
	{
		vec2 const u = {rand_uniform(), rand_uniform()};
		vec3 const p = {center.x + radius * (2 * u.x - 1), center.y + radius * (2 * u.y - 1), 0};
		spawn_particle(p, fluid_type);
	}
}

void scene_structure::spawn_particles_in_sphere(vec3 const &center, float radius, int N, int fluid_type)
{
	for (int k = 0; k < N; ++k)
	{
		vec3 const u = {rand_uniform(), rand_uniform(), rand_uniform()};
		vec3 const p = {center.x + radius * (2 * u.x - 1), center.y + radius * (2 * u.y - 1), center.z + radius * (2 * u.z - 1)};
		spawn_particle(p, fluid_type);
	}
}

void scene_structure::initialize_fluid_classes()
{
	auto water = std::make_shared<fluid_class>("Water", 1.0, 0.0005, vec3{0.1, 0.3, 1.0});
	auto milk = std::make_shared<fluid_class>("Milk", 1.0, 0.0005, vec3{1.0, 1.0, 1.0});
	auto oil = std::make_shared<fluid_class>("Oil", 0.1, 0.05, vec3{1.0, 1.0, 0.3});

	water->soluble_classes.insert(milk);
	milk->soluble_classes.insert(water);

	fluid_classes.push_back(water);
	fluid_classes.push_back(milk);
	fluid_classes.push_back(oil);
}

void scene_structure::initialize_sph()
{
	// Initial particle spacing (relative to h)
	float const c = 1.0f;
	float const h = sph_parameters.h;

	// Fill a square with particles
	particles.clear();
	for (float x = -1.0f + h; x < 1.0f - h; x = x + c * h)
	{
		for (float y = -1.0f + h; y < 1.0f - h; y = y + c * h)
		{
			vec3 pos = {0, 0, 0};
			if (dimension == DIM_2D) {
				pos = {x + h / 8.0 * rand_uniform(), y + h / 8.0 * rand_uniform(), 0};
				spawn_random_type_particle(pos);
			}
			else {
				for (float z = -1.0f + h; z < 1.0f - h; z = z + c * h)
				{
					pos = {x + h / 8.0 * rand_uniform(), y + h / 8.0 * rand_uniform(), z + h / 8.0 * rand_uniform()};
					spawn_random_type_particle(pos);
				}
			}
		}
	}

	if (dimension == DIM_2D) {
		// Initialize the spatial grid
		float const cell_size = 2.0f * h;
		
		int const grid_width = int(ceil(2.0f / cell_size));
		int const grid_height = int(ceil(2.0f / cell_size));

		vec2 const min = {-1.0f, -1.0f};

		grid = spatial_grid(cell_size, grid_width, grid_height, min);

		std::cout << "Grid initialized" << std::endl;
	}
}

void scene_structure::display_frame()
{
	// Set the light to the current position of the camera
	// environment.light = camera_control.camera_model.position();

	if (timer.is_running())
	{
		timer.update(); // update the timer to the current elapsed time
		float const dt = 0.005f * timer.scale;

		if (dimension == DIM_2D) {
			// Update the spatial grid
			grid.clear_grid();
			for (size_t k = 0; k < particles.size(); ++k)
			{
				grid.insert_particle(&particles[k]);
			}

			simulate_2d(dt, particles, grid, sph_parameters);
		}
		else
			simulate_3d(dt, particles, sph_parameters);
	}

	if (gui.display_particles)
	{
		for (int k = 0; k < particles.size(); ++k)
		{
			vec3 const &p = particles[k].p;
			sphere_particle.model.translation = p;
			sphere_particle.material.color = particles[k].color;
			draw(sphere_particle, environment);
		}
	}

	if (gui.display_radius)
	{
		curve_visual.model.scaling = sph_parameters.h;
		for (int k = 0; k < particles.size(); k += 10)
		{
			curve_visual.model.translation = particles[k].p;
			draw(curve_visual, environment);
		}
	}

	if (gui.display_color)
	{
		update_field_color();
		field_quad.texture.update(field);
		draw(field_quad, environment);
	}
}

vec3 color_interpolation(vec3 const &c0, vec3 const &c1, float x0, float x1, float x)
{
	float const alpha = (x - x0) / (x1 - x0);
	return (1 - alpha) * c0 + alpha * c1;
}

vec3 color_clamp(vec3 const &v, float vmin, float vmax)
{
	return {clamp(v.x, vmin, vmax), clamp(v.y, vmin, vmax), clamp(v.z, vmin, vmax)};
}

vec3 scene_structure::get_particle_color(particle_element const &particle)
{
	vec3 res = {1, 1, 1};

	switch (gui.color_type)
	{
	case FLUID_COLOR:
	{
		res = particle.color;
		break;
	}
	case VELOCITY:
	{
		vec3 const v = particle.v;
		res = color_clamp(
			color_interpolation(
				gui.color_min, gui.color_max, gui.threshold_min, gui.threshold_max, norm(v)),
			0.0, 1.0);
		break;
	}
	}
	return res;
}

void scene_structure::update_field_closest(int Nf)
{
	// Iterate over each grid cell in the field
	#pragma omp parallel for
	for (int kx = 0; kx < Nf; ++kx)
	{
		for (int ky = 0; ky < Nf; ++ky)
		{
			vec3 full_color = {0, 0, 0};  // Store the resulting color
			vec3 const p0 = {2.0f * (kx / (Nf - 1.0f) - 0.5f), 2.0f * (ky / (Nf - 1.0f) - 0.5f), 0.0f};
			int cell_index = grid.compute_cell_index(p0);
			float min_dist = sph_parameters.h;

			std::vector<particle_element*> neighbors = grid.get_neighbors(cell_index);

			for (particle_element* particle : neighbors)
			{
				vec3 const& color = get_particle_color(*particle);
				vec3 const& pi = particle->p;
				float const r = norm(pi - p0);

				if (r < min_dist)
				{
					min_dist = r;
					full_color = color;
				}
			}

			field(kx, Nf - 1 - ky) = full_color;
		}
	}
}

void scene_structure::update_field_mean(int Nf)
{
	#pragma omp parallel for
	for (int kx = 0; kx < Nf; ++kx)
	{
		for (int ky = 0; ky < Nf; ++ky)
		{
			vec3 full_color = {0, 0, 0};
			int nb = 0;
			vec3 const p0 = {2.0f * (kx / (Nf - 1.0f) - 0.5f), 2.0f * (ky / (Nf - 1.0f) - 0.5f), 0.0f};

			for (size_t k = 0; k < particles.size(); ++k)
			{
				vec3 const &color = get_particle_color(particles[k]);
				vec3 const &pi = particles[k].p;
				float const r = norm(pi - p0);

				if (r < gui.color_smoothing_radius) {
					full_color += color;
					nb += 1;
				}
			}

			if (nb > 0)
				full_color /= nb;
			field(kx, Nf - 1 - ky) = full_color;
		}
	}
}

void scene_structure::update_field_color()
{
	field.fill({1, 1, 1});
	int const Nf = int(field.dimension.x);

	if (gui.color_type == FLUID_COLOR)
	{
		update_field_closest(Nf);
	}
	else if (gui.color_type == VELOCITY)
	{
		update_field_mean(Nf);
	}
}

void scene_structure::delete_particles_in_disk(vec3 const &center, float radius)
{
	std::vector<particle_element> new_particles;

	for (size_t k = 0; k < particles.size(); ++k)
	{
		vec3 const &p = particles[k].p;
		if (norm(p - center) > radius)
			new_particles.push_back(particles[k]);
	}
	particles = new_particles;
}

void scene_structure::add_radial_force(vec3 const &center, float radius)
{
	#pragma omp parallel
	for (size_t k = 0; k < particles.size(); ++k)
	{
		vec3 const &p = particles[k].p;
		if (norm(p - center) < radius) {
			vec3 force = 100.0f * gui.force_strength * (p - center);
			particles[k].external_forces += force;
		}
	}
}

void scene_structure::mouse_move_event()
{
	if (inputs.mouse.click.left) {
		camera_control.action_mouse_move(environment.camera_view, dimension, cam_projection);
	}
}

void scene_structure::mouse_click_event()
{
	if (inputs.mouse.click.left)
	{ // Movements
		camera_control.action_mouse_move(environment.camera_view, dimension, cam_projection);
	}

	if (inputs.mouse.click.right)
	{ // Special action
		vec2 const cursor = inputs.mouse.position.current;
		vec3 const p = {cursor.x, cursor.y, 0};

		if (gui.right_click_action == SPAWN_PARTICLES)
		{
			if (dimension == DIM_2D){
				spawn_particles_in_disk(p, gui.spawn_particle_radius, gui.spawn_particle_number, gui.spawn_particle_type);
			}
			else {
				// TODO
			}
		}
		else if (gui.right_click_action == REMOVE_PARTICLES)
		{
			if (dimension == DIM_2D) {
				delete_particles_in_disk(p, gui.spawn_particle_radius);
			}
		}
		else if (gui.right_click_action == ADD_FORCE)
		{
			if (dimension == DIM_2D) {
				add_radial_force(p, gui.spawn_particle_radius);
			}
		}
	}
}

void scene_structure::keyboard_event()
{
	camera_control.action_keyboard(environment.camera_view);
}
void scene_structure::idle_frame()
{
	camera_control.idle_frame(environment.camera_view);
}
