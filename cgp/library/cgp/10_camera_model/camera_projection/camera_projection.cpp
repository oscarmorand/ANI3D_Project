#include "camera_projection.hpp"

#include "cgp/01_base/base.hpp"


namespace cgp {

	camera_projection camera_projection::build_perspective(float field_of_view, float aspect_ratio, float depth_min, float depth_max)
	{
		camera_projection res;
		res.type = camera_perspective_type::perspective;
		res.perspective_data.field_of_view = field_of_view;
		res.perspective_data.aspect_ratio = aspect_ratio;
		res.perspective_data.depth_min = depth_min;
		res.perspective_data.depth_max = depth_max;
		return res;
	}
		
	camera_projection camera_projection::build_orthographic(float left, float right, float bottom, float top, float back, float front) 
	{
		camera_projection res;
		res.type = camera_perspective_type::orthographic;
		res.orthographic_data.left = left;
		res.orthographic_data.right = right;
		res.orthographic_data.bottom = bottom;
		res.orthographic_data.top = top;
		res.orthographic_data.z_near = back;
		res.orthographic_data.z_far = front;
		return res;
	}

	void camera_projection::update_aspect_ratio(float aspect_ratio)
	{
		this->aspect_ratio = aspect_ratio;
	}

	mat4 camera_projection::matrix() const
	{
		if (type == camera_perspective_type::perspective)
			return projection_perspective(
				perspective_data.field_of_view,
				perspective_data.aspect_ratio, 
				perspective_data.depth_min, 
				perspective_data.depth_max);
		else
			return projection_orthographic(
				orthographic_data.left*aspect_ratio,
				orthographic_data.right*aspect_ratio, 
				orthographic_data.bottom, 
				orthographic_data.top, 
				orthographic_data.z_near, 
				orthographic_data.z_far);
	}

	mat4 camera_projection::matrix_inverse() const
	{
		if (type == camera_perspective_type::perspective)
			return projection_perspective_inverse(
				perspective_data.field_of_view,
				perspective_data.aspect_ratio, 
				perspective_data.depth_min, 
				perspective_data.depth_max);
		else
			return projection_orthographic_inverse(
				orthographic_data.left*aspect_ratio,
				orthographic_data.right*aspect_ratio, 
				orthographic_data.bottom, 
				orthographic_data.top, 
				orthographic_data.z_near, 
				orthographic_data.z_far);
	}

}