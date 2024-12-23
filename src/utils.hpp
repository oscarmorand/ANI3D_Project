#pragma once

enum color_type_enum {
	FLUID_COLOR,
	VELOCITY,
	DENSITY
};

enum dimension_enum {
	DIM_2D,
	DIM_3D
};

enum right_click_action_enum {
	SPAWN_PARTICLES,
	REMOVE_PARTICLES,
	ADD_RADIAL_FORCE,
	ADD_VORTEX_FORCE,
	ADD_GRAVITY_FORCE
};

enum vortex_direction_enum {
	CLOCKWISE,
	COUNTER_CLOCKWISE
};

enum display_mode_3d_enum {
	MESH,
	POINTS
};