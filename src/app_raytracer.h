#include "def.h"
#include "util.h"
#include "util_optics.h"
#include "material.h"
#include "light.h"
#include "intersect.h"
#include "cornell_box.h"

// ----------------------------------------------------------------------------
// Raytracer
// ----------------------------------------------------------------------------

vec3 background(_in(ray_t) ray)
{
	return vec3(0, 0, 0);
}

void setup_scene()
{
	materials[mat_debug].base_color = vec3(1., 1., 1.);
	materials[mat_debug].metallic = 0.;
	materials[mat_debug].roughness = 0.;
	materials[mat_debug].ior = 1.;
	materials[mat_debug].reflectivity = 0.;
	materials[mat_debug].translucency = 0.;

	setup_cornell_box();

#if 1
	float _sin = sin(u_time);
	float _cos = cos(u_time);
	cb_spheres[cb_sphere_left].origin += vec3(0, abs(_sin), _cos + 1.);
	cb_spheres[cb_sphere_right].origin.z = 0.;// += vec3(0, abs(_cos), _cos);
	lights[0].L.z = 1.5;
#endif
}

void setup_camera(_inout(vec3) eye, _inout(vec3) look_at)
{
	vec2 mouse = u_mouse.x < BIAS ? vec2(0, 0) : 2. * (u_res.xy / u_mouse.xy) - 1.;
	mat3 rot_y = rotate_around_y(mouse.x * 30.);
	eye = mul(rot_y, vec3(0, cb_plane_dist, 2.333 * cb_plane_dist));
	look_at = vec3(0, cb_plane_dist, 0);
}

vec3 illuminate(_in(vec3) eye, _in(hit_t) hit) // TODO: find a way to account for more light types
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
	int i;

	for (i = 0; i < num_cb_planes; ++i) {
		intersect_plane(ray, cb_planes[i], hit);
	}

	for (i = 0; i < num_cb_spheres; ++i) {
		if (cb_spheres[i].material != mat_to_ignore) {
			intersect_sphere(ray, cb_spheres[i], hit);
		}
	}

	return hit;
}

vec3 render(_in(ray_t) primary_ray)
{
	vec3 color = vec3(0, 0, 0);
	vec3 accum = vec3(1, 1, 1);
	ray_t ray = primary_ray;

	for (int i = 0; i < 2; i++) {
		hit_t hit = raytrace_iteration(ray, mat_invalid);

		if (hit.t >= max_dist) {
			color += accum * background(ray);
			break;
		}

		float f = fresnel_factor(1., 1., dot(hit.normal, -ray.direction));
		color += (1. - f) * accum * illuminate(primary_ray.origin, hit);

#if 1 // shadow ray
		if (i == 0) {
			vec3 shadow_line = lights[0].L - hit.origin; // TODO: more light types
			vec3 shadow_dir = normalize(shadow_line);

			ray_t shadow_trace = _begin(ray_t)
				hit.origin + shadow_dir * BIAS,
				shadow_dir
				_end;
			hit_t shadow_hit = raytrace_iteration(shadow_trace, mat_debug);

			if (shadow_hit.t < length(shadow_line)) {
				color *= 0.1;
			}
		}
#endif

		material_t mat = get_material(hit.material_id);
		if (mat.reflectivity > 0.) {
			accum *= f;
			vec3 reflect_dir = normalize(reflect(hit.normal, ray.direction));
			ray.origin = hit.origin + reflect_dir * BIAS;
			ray.direction = reflect_dir;
		} else {
			break;
		}
	}

	return color;
}

#include "main.h"