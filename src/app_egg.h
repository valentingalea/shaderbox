#include "def.h"
#include "util.h"
#include "IK.h"
#include "sdf.h"

// ----------------------------------------------------------------------------
// Vectorpark Egg
// ----------------------------------------------------------------------------

vec3 background(_in(ray_t) ray)
{
	return vec3(.1, .1, .7);
}

void setup_scene()
{
#define mat_debug 0
#define mat_egg 1
#define mat_bike 2
#define mat_ground 3
}

void setup_camera(_inout(vec3) eye, _inout(vec3) look_at)
{
	eye = vec3(.0, .25, 5.25);
	look_at = vec3(.0, .25, .0);
}

vec3 illuminate(_in(hit_t) hit)
{
	if (hit.material_id == mat_ground) return vec3(13. / 255., 104. / 255., 0. / 255.);
	if (hit.material_id == mat_egg) return vec3(0.9, 0.95, 0.95);
	if (hit.material_id == mat_bike) return vec3(.2, .2, .2);
	return vec3(1, 1, 1);
}

#define BEZIER
vec2 sdf(_in(vec3) P)
{
	vec3 p = mul(rotate_around_y(u_time * -100.0), P)
		- vec3(0, 0.5, 3.5);

	int material = mat_egg;

	float egg_y = 0.65;
#if 1
	float egg_m = sd_sphere(p - vec3(0, egg_y, 0), 0.475);
	float egg_b = sd_sphere(p - vec3(0, egg_y - 0.45, 0), 0.25);
	float egg_t = sd_sphere(p - vec3(0, egg_y + 0.45, 0), 0.25);
	float egg_1 = op_blend(egg_m, egg_b, .5);
	float egg_2 = op_blend(egg_1, egg_t, .5);
	vec2 egg = vec2(egg_2, material);
#else
	float s = 1.55;
	mat3 scale = mat3(
		s, 0, 0,
		0, 1, 0,
		0, 0, 1);
	mat3 iscale = mat3(
		1./s, 0, 0,
		0, 1./s, 0,
		0, 0, 1.);
	vec2 egg = vec2(
		sd_sphere(iscale * (scale * (p - vec3(0, egg_y, 0))), 0.475),
		material);
#endif

	vec3 wheel_pos = vec3(0, 1.2, 0);
	float pedal_radius = 0.3;
	float pedal_speed = 400.;
	float pedal_off = 0.2;

	mat3 rot_z = rotate_around_z(-u_time * pedal_speed);
	vec3 left_foot_pos = wheel_pos + mul(rot_z, vec3(0., pedal_radius, pedal_off));

	rot_z = rotate_around_z(-u_time * pedal_speed);
	vec3 right_foot_pos = wheel_pos + mul(rot_z, vec3(0., -pedal_radius, -pedal_off));

	vec3 side = vec3(0, 0, pedal_off);
	float femur = 0.8;
	float tibia = 0.75;
	float thick = .05;

	vec3 pelvis = vec3(0, 0., 0) + side;
	vec3 knee_l = ik_solver(pelvis, left_foot_pos, femur, tibia);
#ifndef BEZIER
	vec2 left_leg_a = vec2(
		sd_cylinder(p + pelvis, vec3(0., 0., 0.), knee_l - side, thick),
		material);
	vec2 left_leg_b = vec2(
		sd_cylinder(p + knee_l, vec3(0., 0., 0.), left_foot_pos - knee_l, thick),
		material);
#endif

	pelvis = vec3(0, 0., 0) - side;
	vec3 knee_r = ik_solver(pelvis, right_foot_pos, femur, tibia);
#ifndef BEZIER
	vec2 right_leg_a = vec2(
		sd_cylinder(p + pelvis, vec3(0., 0., 0.), knee_r + side, thick),
		material);
	vec2 right_leg_b = vec2(
		sd_cylinder(p + knee_r, vec3(0., 0., 0.), right_foot_pos - knee_r, thick),
		material);
#endif

	vec2 legs = op_add(
#ifndef BEZIER
		vec2(op_blend(left_leg_a.x, left_leg_b.x, .01), material),
		op_add(right_leg_a, right_leg_b)
#else
		vec2(
		sd_bezier(-(vec3(0., 0., 0.) + side), -knee_l, -left_foot_pos, p, thick).x,
		material),
		vec2(
		sd_bezier(-(vec3(0., 0., 0.) - side), -knee_r, -right_foot_pos, p, thick).x,
		material)
#endif
	);

	vec3 left_toe = normalize(vec3(left_foot_pos.y - knee_l.y, knee_l.x - left_foot_pos.x, 0));
	vec2 left_foot = vec2(
		sd_cylinder(p + left_foot_pos, vec3(0., 0., 0.), left_toe / 8., thick),
		material);

	vec3 right_toe = normalize(vec3(right_foot_pos.y - knee_r.y, knee_r.x - right_foot_pos.x, 0));
	vec2 right_foot = vec2(
		sd_cylinder(p + right_foot_pos, vec3(0., 0., 0.), right_toe / 8., thick),
		material);

	vec2 feet = op_add(left_foot, right_foot);

	vec2 bike = vec2(
		sd_torus(p + wheel_pos, 1., .03),
		mat_bike);

	vec2 ground = vec2(
		sd_plane(P, vec3(0., 1., 0.), wheel_pos.y + 0.5),
		mat_ground);

	vec2 _1 = op_add(feet, bike);
	vec2 _2 = op_add(egg, _1);
	vec2 _3 = op_add(legs, _2);
	return op_add(ground, _3);
}

vec3 sdf_normal(_in(vec3) p)
{
	float dt = 0.05;
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

float shadowmarch(_in(ray_t) ray)
{
	const int steps = 20;
	const float end = 10.;
	const float penumbra_factor = 15.;
	const float darkest = 0.1;

	float t = 0.;
	float umbra = 1.;
	for (int i = 0; i < steps; i++) {
		vec3 p = ray.origin + ray.direction * t;
		vec2 d = sdf(p);

		if (t > end) break;
		if (d.x < EPSILON) {
			return darkest;
		}

		t += d.x;
		
		// from http://iquilezles.org/www/articles/rmshadows/rmshadows.htm
		umbra = min(umbra, penumbra_factor * d.x / t);
	}

	return umbra;
}

_mutable(float) depth = -max_dist;

vec3 render(_in(ray_t) ray)
{
	const int steps = 80;
	const float end = 15.;

	float t = 0.;
	for (int i = 0; i < steps; i++) {
		vec3 p = ray.origin + ray.direction * t;
		vec2 d = sdf(p);

		if (t > end) break;
		if (d.x < EPSILON) {
			hit_t h = _begin(hit_t)
				t, // ray length at impact
				int(d.y), // material id
				vec3(0, 0, 0), // sdf_normal(p),
				p // point of impact				
			_end;

			if (h.material_id == mat_egg || h.material_id == mat_bike) {
				depth = max(depth, p.z);
			}

			float s = 1.;
#if 1 // soft shadows
			if (int(d.y) == mat_ground) {
				vec3 sh_dir = vec3(0, 1, 1);
				ray_t sh_ray = _begin(ray_t)
					p + sh_dir * 0.05, sh_dir
				_end;
				s = shadowmarch(sh_ray);
			}
#endif

			return illuminate(h) * s;
		}

		t += d.x;
	}

	return background(ray);
}

void mainImage(
	_out(vec4) fragColor,
#ifdef SHADERTOY
	vec2 fragCoord
#else
	_in(vec2) fragCoord
#endif
){
	vec2 aspect_ratio = vec2(u_res.x / u_res.y, 1);
	float fov = tan(radians(30.0));

	vec3 final_color = vec3(0, 0, 0);

	vec3 eye, look_at;
	setup_camera(eye, look_at);

	setup_scene();

	vec2 point_ndc = fragCoord.xy / u_res.xy;
#ifdef HLSL
		point_ndc.y = 1. - point_ndc.y;
#endif
	vec3 point_cam = vec3(
		(2.0 * point_ndc - 1.0) * aspect_ratio,// * fov,
		-1.0);
	ray_t ray = get_primary_ray(point_cam, eye, look_at);

	final_color += render(ray);

#if 1
	// from https://www.shadertoy.com/view/4sjGzc
#define BAR_SEPARATION 0.6
#define BAR_WIDTH 0.05
#define BAR_DEPTH 1.
#define BAR_COLOR vec3(.6, .6, .6)
	float bar_factor = 1.0 - smoothstep(0.0, 0.01, abs((abs(point_cam.x) - BAR_SEPARATION)) - BAR_WIDTH);
	float depth_factor = 1. - step(BAR_DEPTH, depth);
	final_color = mix(final_color, BAR_COLOR, bar_factor * depth_factor);
#endif

	fragColor = vec4(linear_to_srgb(final_color), 1.);
}
