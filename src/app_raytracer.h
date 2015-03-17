#include "def.h"
#include "util.h"
#include "util_optics.h"
#include "material.h"
#include "light.h"
#include "intersect.h"

#define num_planes 6
plane_t planes[num_planes];

#define num_spheres 3
sphere_t spheres[num_spheres];

vec3 background(_in(ray_t) ray)
{
	return vec3(0);
}

#include "cornell_box.h"
void setup_scene()
{
	materials[mat_debug] = material_t _begin vec3(1., 1., 1.), 0., 0., 1., 0., 0. _end;

	setup_cornell_box();
}

void setup_camera(_inout(vec3) eye, _inout(vec3) look_at)
{
	vec2 mouse = iMouse.x < BIAS ? vec2(0) : 2. * (iResolution.xy / iMouse.xy) - 1.;
	mat3 rot_y = rotate_around_y(mouse.x * 30.);
	eye = rot_y * vec3(0, cb_plane_dist, 2.333 * cb_plane_dist);
	look_at = vec3(0, cb_plane_dist, 0);
}

vec3 illuminate(_in(hit_t) hit) // TODO: find a way to account for more light types
{
	material_t mat = get_material(hit.material_id);

	// special case for debug stuff - just solid paint it
	if (hit.material_id == mat_debug) {
		return materials[mat_debug].base_color;
	}

	vec3 accum = ambient_light; // really cheap equivalent for indirect light

	vec3 V = normalize(eye - hit.origin); // view direction
	vec3 L = get_light_direction(lights[0], hit);

	// TODO: more lights
#if 0
		accum += illum_blinn_phong(V, L, hit, mat);
#else
		accum += illum_cook_torrance(V, L, hit, mat);
#endif

	return accum;
}

hit_t raytrace_iteration(_in(ray_t) ray, _in(int) mat_to_ignore)
{
	hit_t hit = no_hit;

	for (int i = 0; i < num_planes; ++i) {
		intersect_plane(ray, planes[i], hit);
	}

	for (int i = 0; i < num_spheres; ++i) {
		if (spheres[i].material != mat_to_ignore) {
			intersect_sphere(ray, spheres[i], hit);
		}
	}

	return hit;
}

vec3 render(_in(ray_t) ray)
{
	hit_t hit = raytrace_iteration(ray, mat_invalid);

	if (hit.t >= max_dist) {
		return background(ray);
	}

	vec3 color = illuminate(hit);

#if 1 // shadow ray
	vec3 sh_line = lights[0].L - hit.origin; // TODO: more light types
	vec3 sh_dir = normalize(sh_line);
	ray_t sh_trace = ray_t _begin
		hit.origin + sh_dir * BIAS,
		sh_dir
	_end;
	hit_t sh_hit = raytrace_iteration(sh_trace, mat_debug);
	if (sh_hit.t < length(sh_line)) {
		color *= 0.1;
	}
#endif

	return color;
}

#include "main.h"