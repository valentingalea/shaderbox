// ----------------------------------------------------------------------------
// Cornell Box definition
// http://www.graphics.cornell.edu/online/box/
// ----------------------------------------------------------------------------

// reuses materials array

#define num_cb_planes 6
_mutable(plane_t) cb_planes[num_cb_planes];

#define num_cb_spheres 3
_mutable(sphere_t) cb_spheres[num_cb_spheres];

void setup_material(
	_inout(material_t) mat,
	_in(vec3) diffuse,
	_in(float) metallic,
	_in(float) roughness
){
	mat.base_color = diffuse;
	mat.metallic = metallic;
	mat.roughness = roughness;
	mat.ior = 1.;
	mat.reflectivity = 0.;
	mat.translucency = 0.;
}

void setup_plane(
	_inout(plane_t) p,
	_in(vec3) n,
	_in(float) d,
	_in(int) mat_id
){
	p.direction = n;
	p.distance = d;
	p.material = mat_id;
}

void setup_cornell_box()
{
#define cb_mat_white 1
#define cb_mat_red 2
#define cb_mat_blue 3
#define cb_mat_reflect 4
#define cb_mat_refract 5
#define cb_mat_green 6
	setup_material(materials[cb_mat_white], vec3(0.7913, 0.7913, 0.7913), .0, .5);
	setup_material(materials[cb_mat_red], vec3(0.6795, 0.0612, 0.0529), 0., .5);
	setup_material(materials[cb_mat_blue], vec3(0.1878, 0.1274, 0.4287), 0., .5);
	setup_material(materials[cb_mat_reflect], vec3(0.95, 0.64, 0.54), 1., .1);
	materials[cb_mat_reflect].reflectivity = 1.;
	setup_material(materials[cb_mat_refract], vec3(1., 0.77, 0.345), 1., .05);
	materials[cb_mat_refract].reflectivity = 1.;
	materials[cb_mat_refract].translucency = 0.;
	materials[cb_mat_refract].ior = 1.333;

#define cb_plane_ground 0
#define cb_plane_behind 1
#define cb_plane_front 2
#define cb_plane_ceiling 3
#define cb_plane_left 4
#define cb_plane_right 5
#define cb_plane_dist 2.
	setup_plane(cb_planes[cb_plane_ground], vec3(0, -1, 0), 0., cb_mat_white);
	setup_plane(cb_planes[cb_plane_ceiling], vec3(0, 1, 0), 2. * cb_plane_dist, cb_mat_white);
	setup_plane(cb_planes[cb_plane_behind], vec3(0, 0, -1), -cb_plane_dist, cb_mat_white);
	setup_plane(cb_planes[cb_plane_front], vec3(0, 0, 1), cb_plane_dist, cb_mat_white);
	setup_plane(cb_planes[cb_plane_left], vec3(1, 0, 0), cb_plane_dist, cb_mat_red);
	setup_plane(cb_planes[cb_plane_right], vec3(-1, 0, 0), -cb_plane_dist, cb_mat_blue);

#define cb_sphere_light 0
#define cb_sphere_left 1
#define cb_sphere_right 2
	cb_spheres[cb_sphere_light].origin = vec3(0, 2.5 * cb_plane_dist + 0.4, 0);
	cb_spheres[cb_sphere_light].radius = 1.5;
	cb_spheres[cb_sphere_light].material = mat_debug;
	cb_spheres[cb_sphere_left].origin = vec3(0.75, 1, -0.75);
	cb_spheres[cb_sphere_left].radius = 0.75;
	cb_spheres[cb_sphere_left].material = cb_mat_reflect;
	cb_spheres[cb_sphere_right].origin = vec3(-0.75, 0.75, 0.75);
	cb_spheres[cb_sphere_right].radius = 0.75;
	cb_spheres[cb_sphere_right].material = cb_mat_refract;

	lights[0].type = LIGHT_POINT;
	lights[0].L = vec3(0, 2. * cb_plane_dist - 0.2, 0);
	lights[0].color = vec3(1., 1., 1.);
}

