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
	
	image_structure forest_image = image_load_file(project::path + "assets/images/skybox_01.jpg");
	std::vector<image_structure> image_grid = image_split_grid(forest_image, 4, 3);
	skybox.initialize_data_on_gpu();
	skybox.texture.initialize_cubemap_on_gpu(image_grid[1], image_grid[7], image_grid[5], image_grid[3], image_grid[10], image_grid[4]);

	shader_environment_map.load(project::path + "shaders/mesh/environment_map.vert.glsl", project::path + "shaders/mesh/environment_map.frag.glsl");
	environment.uniform_generic.uniform_mat3["skybox_rotation"] = mat3::build_identity();
	environment.default_expected_uniform = false;
	vec3 const length = {2.4,2.4,2.4};
	implicit_surface.set_domain(50, length);
}

void scene_structure::spawn_particle(vec3 const &pos, int fluid_type, vec3 const &velocity)
{
	particle_element particle;
	particle.p = pos;
	particle.v = velocity;
	particle.f = {0, 0, 0};

	auto const& fluid_class = fluid_classes[fluid_type];
	particle.m = fluid_class->base_m;
	particle.nu = fluid_class->base_nu;
	particle.color = fluid_class->base_color;

	particle.fluid_type = fluid_class;

	particles.push_back(particle);
	gui.nb_particles += 1;
}

void scene_structure::spawn_particle(vec3 const &pos, int fluid_type)
{
	spawn_particle(pos, fluid_type, {0, 0, 0});
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
		vec3 const velocity = p - u;
		spawn_particle(p, fluid_type, velocity * 3);
	}
}

void scene_structure::initialize_fluid_classes()
{
	auto water = std::make_shared<fluid_class>("Water", 1.0, 0.0005, vec3{0.1, 0.3, 1.0});
	auto milk = std::make_shared<fluid_class>("Milk", 1.0, 0.0005, vec3{1.0, 1.0, 1.0});
	auto oil = std::make_shared<fluid_class>("Oil", 0.92, 0.05, vec3{1.0, 1.0, 0.3});

	water->soluble_classes.insert(milk);
	milk->soluble_classes.insert(water);

	fluid_classes.push_back(water);
	fluid_classes.push_back(milk);
	fluid_classes.push_back(oil);
}

void scene_structure::initialize_sph()
{
	// Initial particle spacing (relative to h)
	float const c = 3.0f;
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

	float const cell_size = 2.0f * h;
		
	int const grid_width = int(ceil(2.0f / cell_size));
	int const grid_height = int(ceil(2.0f / cell_size));

	if (dimension == DIM_2D) {
		vec2 const min = {-1.0f, -1.0f};

		grid_2d = spatial_grid_2d(cell_size, grid_width, grid_height, min);
	}
	else {
		int const grid_depth = int(ceil(2.0f / cell_size));

		vec3 const min = {-1.000001f, -1.000001f, -1.000001f};

		grid_3d = spatial_grid_3d(cell_size, grid_width, grid_height, grid_depth, min);
	}
	std::cout << "Grid initialized" << std::endl;
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
			grid_2d.clear_grid();
			for (size_t k = 0; k < particles.size(); ++k)
			{
				grid_2d.insert_particle(&particles[k]);
			}

			simulate_2d(dt, particles, grid_2d, sph_parameters);
		}
		else {
			grid_3d.clear_grid();
			for (size_t k = 0; k < particles.size(); ++k)
			{
				grid_3d.insert_particle(&particles[k]);
			}
			
			simulate_3d(dt, particles, grid_3d, sph_parameters);
		}

		if (inputs.mouse.click.right)
		{ // Special action
			right_click();
		}
	}

	if (gui.display_particles)
	{	
		//Create a field that represents each particle force.
		for (int k = 0; k < particles.size(); ++k)
		{
			vec3 const &p = particles[k].p;
			sphere_particle.model.translation = p;
			sphere_particle.model.scaling = gui.particle_radius_ratio * sph_parameters.h;
			sphere_particle.material.color = get_particle_color(particles[k]);
			draw(sphere_particle, environment);
		}
	}
	if (gui.display_skybox)
	{
		if (dimension == DIM_3D) {
			glDepthMask(GL_FALSE);
			draw(skybox, environment);
			glDepthMask(GL_TRUE);
		}
	}
	if (gui.display_mesh)
	{
		if (dimension == DIM_3D) {
			glDisable(GL_DEPTH_TEST);
			implicit_surface.update_field(field_function, grid_3d, gui.alpha, gui.influence_radius_MC, gui.isovalue_MC, gui.reflectivness);
			implicit_surface.drawable_param.shape.shader = shader_environment_map;
			implicit_surface.drawable_param.shape.supplementary_texture["image_skybox"] = skybox.texture;
			draw(implicit_surface.drawable_param.shape, environment);
			glEnable(GL_DEPTH_TEST);
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
	gui.nb_particles = particles.size();
}

void scene_structure::add_vortex_force(vec3 const &center, float radius, float strength)
{
	#pragma omp parallel
	for (size_t k = 0; k < particles.size(); ++k)
	{
		vec3 const &p = particles[k].p;
		vec3 const diff = center - p;
		if (norm(diff) < radius) {
			vec3 dir = cross(vec3(0,0,-1), normalize(diff));
			dir = dir * (gui.vortex_direction == CLOCKWISE ? -1.0f : 1.0f);
			vec3 force = strength * dir;
			particles[k].external_forces += force;
		}
	}
}

void scene_structure::add_radial_force(vec3 const &center, float radius, float strength)
{
	#pragma omp parallel
	for (size_t k = 0; k < particles.size(); ++k)
	{
		vec3 const &p = particles[k].p;
		vec3 diff = p - center;
		if (norm(diff) < radius) {
			diff = normalize(diff);
			vec3 force = strength * diff;
			particles[k].external_forces += force;
		}
	}
}

void scene_structure::add_gravity_force(vec3 const &center, float radius, float strength)
{
	#pragma omp parallel
	for (size_t k = 0; k < particles.size(); ++k)
	{
		vec3 const &p = particles[k].p;
		vec3 diff = p - center;
		if (norm(diff) < radius) {
			vec3 force = 10.0f * strength * (-diff);
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

base_plan scene_structure::get_most_orthogonal_plan(vec3 const &dir) 
{
	base_plan plan_x = {vec3{1, 0, 0}};
	base_plan plan_y = {vec3{0, 1, 0}};
	base_plan plan_z = {vec3{0, 0, 1}};

	float angle_x = dot(dir, plan_x.normal);
	float angle_y = dot(dir, plan_y.normal);
	float angle_z = dot(dir, plan_z.normal);

	if (angle_x > angle_y && angle_x > angle_z)
		return plan_x;
	else if (angle_y > angle_x && angle_y > angle_z)
		return plan_y;
	else
		return plan_z;
}

void scene_structure::right_click()
{
	vec2 const cursor = inputs.mouse.position.current;
	vec3 const p = {cursor.x, cursor.y, 0};

	if (gui.right_click_action == SPAWN_PARTICLES)
		spawn_particles_in_disk(p, gui.spawn_particle_radius, gui.spawn_particle_number, gui.spawn_particle_type);
	else if (gui.right_click_action == REMOVE_PARTICLES)
		delete_particles_in_disk(p, gui.spawn_particle_radius);
	else if (gui.right_click_action == ADD_RADIAL_FORCE)
		add_radial_force(p, gui.spawn_particle_radius, gui.force_strength);
	else if (gui.right_click_action == ADD_VORTEX_FORCE)
		add_vortex_force(p, gui.spawn_particle_radius, gui.force_strength);
	else if (gui.right_click_action == ADD_GRAVITY_FORCE)
		add_gravity_force(p, gui.spawn_particle_radius, gui.force_strength);
}

void scene_structure::mouse_click_event()
{
	if (inputs.mouse.click.left)
	{ // Movements
		camera_control.action_mouse_move(environment.camera_view, dimension, cam_projection);
	}

	if (inputs.mouse.click.right)
	{ // Special action
		right_click();
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
