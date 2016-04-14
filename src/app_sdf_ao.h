#include "def.h"
#include "util.h"
#include "sdf.h"

// ----------------------------------------------------------------------------
// Distance Fields Ambient Occlusion
// ----------------------------------------------------------------------------

vec3 background(_in(ray_t) ray)
{
	return vec3(.1, .1, .7);
}

#define mat_debug	0
#define mat_ground	1
#define mat_pipe	2
#define mat_bottom	3
#define mat_deck	4
#define mat_coping	5
#define mat_count	6
_mutable(vec3) materials[mat_count];

vec3 get_material(_in(int) index)
{
	vec3 mat;
	for (int i = 0; i < mat_count; ++i) {
		if (i == index) {
			mat = materials[i];
			break;
		}
	}
	return mat;
}

void setup_scene()
{
	materials[mat_debug] = vec3(1, 1, 1);
	materials[mat_ground] = vec3(0, 1, 0);
	materials[mat_pipe] = vec3(.6, .6, .6);
	materials[mat_bottom] = materials[mat_pipe];
	materials[mat_pipe] = materials[mat_pipe];
	materials[mat_coping] = vec3(.9, .9, .9);
}

void setup_camera(_inout(vec3) eye, _inout(vec3) look_at)
{
	mat3 rot = rotate_around_y (u_time * 50.);
	eye = mul(rot, vec3(0, 3, 5));
	look_at = vec3(0, 0, 0);
}

vec3 illuminate(_in(vec3) eye, _in(hit_t) hit)
{
#if 0 // debug: output the raymarching steps
	return vec3(hit.normal);
#endif

	vec3 base_color = get_material(hit.material_id);

	vec3 V = normalize(eye - hit.origin); // view direction
	vec3 L = normalize (vec3 (0, 1, 0)); //get_light_direction(lights[0], hit);
	
	return max(0., dot(L, hit.normal))
	* base_color + vec3(.01, .01, .01);
}

_constant(vec3) size = vec3(1.3, 1., 1.25);

vec2 sdf_pipe(_in(vec3) pos)
{
	// origin
	// ramp(box and cylinder)
	vec3 p = pos - vec3(0, size.y, 0);

	float b = sd_box(p, size);

	p -= vec3(.7, .5, 0);
	p = mul(p, rotate_around_x(-90.));
	float c = sd_y_cylinder(p,
		size.y + .55, // radius
		2. * size.z + .1); // height

	vec2 pipe = vec2(
		op_sub(b, c),
		mat_pipe);

	// revert
	// coping bars
	p = pos - vec3(0, size.y, 0);

	p -= vec3(-size.x + .525, size.y, 0);
	p = mul(p, rotate_around_x(-90.));
	vec2 coping = vec2(
		sd_y_cylinder(p,
			.025, // radius
			2. * size.z), // height
		mat_coping);

	// revert
	// the deck railing
	p = pos - vec3(0, size.y * 2., 0);

	float rail = sd_box(
		p + vec3(size.x, -.25, 0),
		vec3(.025, .05, size.z));

	const vec3 B = vec3(.025, .125, .025);
	const float H = -.125;
	float bar_1 = sd_box(p + vec3(size.x, H, 0), B);
	float bar_2 = sd_box(p + vec3(size.x, H, size.z / 2.), B);
	float bar_3 = sd_box(p + vec3(size.x, H, size.z), B);
	float bar_4 = sd_box(p + vec3(size.x, H, -size.z / 2.), B);
	float bar_5 = sd_box(p + vec3(size.x, H, -size.z), B);
	float b_a = op_add(bar_1, bar_2);
	float b_b = op_add(b_a, bar_3);
	float b_c = op_add(bar_4, bar_5);
	float b_d = op_add(b_b, b_c);
	float bars = b_d;

	vec2 railing = vec2(
		op_add(rail, bars),
		mat_deck);
	vec2 deck = op_add(railing, coping);

	return op_add(pipe, deck);
}

vec2 sdf(_in(vec3) pos)
{
	// NOTE: everything is centered around origin
	// change coord frame by offseting
	// with inverse then doing
	// the opposite before next
	// effectively doing push/pop

	// NOTE: all measurements are in halfs
	// due to the above

	const float B = .15;
	vec3 p = pos - vec3(0, B, 0);

	vec2 bottom = vec2(
		sd_box(p, vec3(2.25 * size.x, B / 2., size.z)),
		mat_bottom);

	vec2 pipe1 = sdf_pipe(p + vec3(1.25 * size.x, 0, 0));

	p -= vec3(1.25 * size.x, 0, 0);
	p = mul(p, rotate_around_y(180.));
	vec2 pipe2 = sdf_pipe(p);

	vec2 pipe = op_add(pipe1, pipe2);

	vec2 ref = vec2(
		sd_box(pos, vec3(.025, 15, .025)),
		mat_debug);

	vec2 ground = vec2(
		sd_plane(pos, vec3(0, 1, 0), 0.),
		mat_ground);

	vec2 g = ground; // op_add(ground, ref);
	vec2 b = op_add(pipe, bottom);
	return op_add(b, g);
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

vec3 sdf_ao(_in(hit_t) hit)
{
	const float dt = .5;
	const int steps = 5;
	
	float d = 0.;
	float occlusion = 0.;
	for (float i = 1.; i <= float(steps); i += 1.) {
		vec3 p = hit.origin + dt * i * hit.normal;
		d = sdf (p).x;
		
		occlusion += 1. / pow(2., i) * (dt * i - d);
	}
	
	float c = 1. - clamp(occlusion, 0., 1.);
	return vec3 (c, c, c);
}

#define EPSILON 0.01

vec3 render(_in(ray_t) ray)
{
	const int steps = 70;
	const float end = 20.;

	float t = 0.;
	for (int i = 0; i < steps; i++) {
		vec3 p = ray.origin + ray.direction * t;
		vec2 d = sdf(p);

		//if (t > end) break;
		if (d.x < .005) {
			hit_t h = _begin(hit_t)
				t, // ray length at impact
				int(d.y), // material id
				sdf_normal(p),
				p // point of impact				
			_end;

			float ambient = sdf_ao (h).x;
			return illuminate(ray.origin, h) * ambient;
		}

		t += d.x;
	}

	return background(ray);
}

#include "main.h"