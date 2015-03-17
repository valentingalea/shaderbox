void setup_cornell_box()
{
#define cb_mat_white 1
#define cb_mat_red 2
#define cb_mat_blue 3
#define cb_mat_reflect 4
#define cb_mat_refract 5
#define cb_mat_green 6
	materials[cb_mat_white] = material_t _begin
		vec3(0.7913), .0, .5, 1., 0., 0.
	_end;
	materials[cb_mat_red] = material_t _begin
		vec3(0.6795, 0.0612, 0.0529),
		0., .5, 1., 0., 0.
	_end;
	materials[cb_mat_blue] = material_t _begin
		vec3(0.1878, 0.1274, 0.4287),
		0., .5, 1., 0., 0.
	_end;
	materials[cb_mat_reflect] = material_t _begin
		vec3(0.95, 0.64, 0.54),
		1., .1, 1.0, 0., 0.
	_end;
	materials[cb_mat_refract] = material_t _begin
		vec3(1., 0.77, 0.345),
		1., .05, 1.333, 0., 1.
	_end;

#define cb_plane_ground 0
#define cb_plane_behind 1
#define cb_plane_front 2
#define cb_plane_ceiling 3
#define cb_plane_left 4
#define cb_plane_right 5
#define cb_plane_dist 2.
	planes[cb_plane_ground] = plane_t _begin vec3(0, -1, 0), 0., cb_mat_white _end;
	planes[cb_plane_ceiling] = plane_t _begin vec3(0, 1, 0), 2. * cb_plane_dist, cb_mat_white _end;
	planes[cb_plane_behind] = plane_t _begin vec3(0, 0, -1), -cb_plane_dist, cb_mat_white _end;
	planes[cb_plane_front] = plane_t _begin vec3(0, 0, 1), cb_plane_dist, cb_mat_white _end;
	planes[cb_plane_left] = plane_t _begin vec3(1, 0, 0), cb_plane_dist, cb_mat_red _end;
	planes[cb_plane_right] = plane_t _begin vec3(-1, 0, 0), -cb_plane_dist, cb_mat_blue _end;

#define cb_sphere_light 0
#define cb_sphere_left 1
#define cb_sphere_right 2
	spheres[cb_sphere_light] = sphere_t _begin vec3(0, 2.5 * cb_plane_dist + 0.4, 0), 1.5, mat_debug _end;
	spheres[cb_sphere_left] = sphere_t _begin vec3(0.75, 1, -0.75), 0.5, cb_mat_reflect _end;
	spheres[cb_sphere_right] = sphere_t _begin vec3(-0.75, 0.75, 0.75), 0.75, cb_mat_refract _end;

	lights[0] = light_t _begin
		LIGHT_POINT,
		vec3(0, 2. * cb_plane_dist - 0.2, 0),
		vec3(1., 1., 1.)
	_end;
}