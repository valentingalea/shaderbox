#include "def.h"
#include "util.h"
#include "sdf.h"
#include "material.h"
#include "util_optics.h"
#include "light.h"

vec3 background(_in(ray_t) ray)
{
	return vec3(.1, .1, .7);
}

void setup_scene()
{
#define mat_debug 0
#define mat_ground 0
#define mat_phong 0
	materials[mat_debug] = material_t _begin vec3(1., 1., 1.), 0., 0., 1., 0., 0. _end;
}

void setup_camera(_inout(vec3) eye, _inout(vec3) look_at)
{
	mat3 rot = rotate_around_y (iGlobalTime * 50.);
	eye = rot * vec3(2, 2, 4);
	look_at = vec3(0);
}

vec3 illuminate(_in(hit_t) hit)
{
#if 0 // debug: output the raymarching steps
	return vec3(hit.normal);//.material_param);
#endif

	material_t mat = get_material(hit.material_id);

	vec3 accum = ambient_light; // really cheap equivalent for indirect light

	vec3 V = normalize(eye - hit.origin); // view direction
	vec3 L = vec3 (1, 1, 0); //get_light_direction(lights[0], hit);

#if 1
		accum += illum_blinn_phong(V, L, hit, mat);
#else
		accum += illum_cook_torrance(V, L, hit, mat);
#endif

	return accum;
}

vec2 sdf(_in(vec3) p)
{
	vec2 box = vec2 (
		sd_box (p + vec3 (1, 0, 0), vec3 (1)),
		mat_phong);
		
	vec2 globe = vec2 (
		sd_sphere (p - vec3 (1, 0, 0), 1.),
		mat_phong);
		
	vec2 ground = vec2 (
		sd_plane (p, vec3 (0, 1, 0), 1),
		mat_ground);

	return op_add (
		op_add (box, globe),
		ground);
}

vec3 sdf_normal(_in(vec3) p)
{
	float dt = 0.001;
	vec3 x = vec3(dt, 0, 0);
	vec3 y = vec3(0, dt, 0);
	vec3 z = vec3(0, 0, dt);
	return normalize(vec3(
		sdf(p + x).r - sdf(p - x).r,
		sdf(p + y).r - sdf(p - y).r,
		sdf(p + z).r - sdf(p - z).r
	));
}

#define EPSILON 0.001

vec3 render(_in(ray_t) ray)
{
	const int steps = 50;
	const float end = 10.;

	float t = 0.;
	for (int i = 0; i < steps; i++) {
		vec3 p = ray.origin + ray.direction * t;
		vec2 d = sdf(p);

		if (t > end) break;
		if (d.x < EPSILON) {
			hit_t h = hit_t _begin
				t, // ray length at impact
				int(d.y), // material id
				float(i) / float(steps), // material custom param
				sdf_normal(p),
				p // point of impact				
			_end;

			return illuminate(h);
		}

		t += d.x;
	}

	return background(ray);
}

#include "main.h"