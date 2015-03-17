#include "def.h"
#include "util.h"
#include "IK.h"
#include "sdf.h"

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
	eye = vec3(1, 2, 5);
	look_at = vec3(0);
}

vec3 illuminate(_in(hit_t) hit)
{
#if 0 // debug: output the raymarching steps
	return vec3(hit.material_param);
#endif

	if (hit.material_id == mat_ground) return vec3(13. / 255., 104. / 255., 0. / 255.);
	if (hit.material_id == mat_egg) return vec3(0.95);
	if (hit.material_id == mat_bike) return vec3(.2);
	return vec3(1);
}

vec2 sdf(_in(vec3) P)
{
	vec3 p = rotate_around_y(iGlobalTime * -50.0) * P -
		vec3(0, 0.5, 1.75);

	int material = mat_egg;

	float egg_y = 0.65;
#if 0
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
	float pedal_speed = 300.;
	float pedal_off = 0.2;

	mat3 rot_z = rotate_around_z(-iGlobalTime * pedal_speed);
	vec3 left_foot_pos = wheel_pos + rot_z * vec3(0, pedal_radius, pedal_off);

	rot_z = rotate_around_z(-iGlobalTime * pedal_speed);
	vec3 right_foot_pos = wheel_pos + rot_z * vec3(0, -pedal_radius, -pedal_off);

	vec3 side = vec3(0, 0, pedal_off);
	float femur = 0.8;
	float tibia = 0.75;
	float thick = .05;

	vec3 pelvis = vec3(0, 0., 0) + side;
	vec3 knee_l = ik_solver(pelvis, left_foot_pos, femur, tibia);
	vec2 left_leg_a = vec2(
		sd_cylinder(p + pelvis, vec3(0), knee_l - side, thick),
		material);
	vec2 left_leg_b = vec2(
		sd_cylinder(p + knee_l, vec3(0), left_foot_pos - knee_l, thick),
		material);

	pelvis = vec3(0, 0., 0) - side;
	vec3 knee_r = ik_solver(pelvis, right_foot_pos, femur, tibia);
	vec2 right_leg_a = vec2(
		sd_cylinder(p + pelvis, vec3(0), knee_r + side, thick),
		material);
	vec2 right_leg_b = vec2(
		sd_cylinder(p + knee_r, vec3(0), right_foot_pos - knee_r, thick),
		material);

	vec2 legs = op_add(
		vec2(op_blend(left_leg_a.x, left_leg_b.x, .01), material),
		op_add(right_leg_a, right_leg_b));

	vec3 left_toe = normalize(vec3(left_foot_pos.y - knee_l.y, knee_l.x - left_foot_pos.x, 0));
	vec2 left_foot = vec2(
		sd_cylinder(p + left_foot_pos, vec3(0), left_toe / 8., thick),
		material);

	vec3 right_toe = normalize(vec3(right_foot_pos.y - knee_r.y, knee_r.x - right_foot_pos.x, 0));
	vec2 right_foot = vec2(
		sd_cylinder(p + right_foot_pos, vec3(0), right_toe / 8., thick),
		material);

	vec2 feet = op_add(left_foot, right_foot);

	vec2 bike = vec2(
		sd_torus(p + wheel_pos, 1., .03),
		mat_bike);

	vec2 ground = vec2(
		sd_plane(P, vec3(0, 1, 0), wheel_pos.y + 0.5),
		mat_ground);

	return op_add(
		ground,
		op_add(legs, op_add(egg, op_add(feet, bike))));
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

vec3 render(_in(ray_t) ray)
{
	const int steps = 75;
	const float end = 15.;

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
				p, // point of impact
				vec3(0) // sdf_normal(p);
			_end;

			float s = 1.;
#if 1 // soft shadows
			if (int(d.y) == mat_ground) {
				vec3 sh_dir = vec3(0, 1, 1);
				ray_t sh_ray = ray_t _begin
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

#include "main.h"
