#include "def.h"
#include "util.h"
#include "intersect.h"

#include "noise_iq.h"
#define noise(x) noise_iq(x)
#include "fbm.h"

_constant(sphere_t) planet = _begin(sphere_t)
	vec3(0, 0, 0), 1, 0
_end;

vec3 background(
	_in(ray_t) eye
){
	_constant(vec3) sun_color = vec3(1., .9, .55);
	float sun_amount = dot(eye.direction, vec3(0, 0, 1));

	vec3 sky = mix(
		vec3(.0, .05, .2),
		vec3(.15, .3, .4),
		1.0 - eye.direction.y);
	sky += sun_color * min(pow(sun_amount, 30.0) * 5.0, 1.0);
	sky += sun_color * min(pow(sun_amount, 10.0) * .6, 1.0);

	return sky;
}

void setup_scene()
{
}

void setup_camera(
	_inout(vec3) eye,
	_inout(vec3) look_at
){
	eye = vec3(0, 0, -2.5);
	look_at = vec3(0, 0, 2);
}

void clouds_map(
	_in(vec3) pos,
	_inout(float) T,
	_inout(vec3) C,
	_inout(float) alpha,
	_in(float) t_step
){
	float dens = fbm(pos * 8, 1.);

	const float coverage = .75;
	const float fuziness = .035;
	//dens *= step(coverage, dens);
	dens *= smoothstep(coverage, coverage + fuziness, dens);

	const float absorbtion = 5.3434;
	float T_i = exp(-absorbtion * dens * t_step);
	T *= T_i;
	C += T // * (exp(hs) / 1.75)
		* dens * t_step;
	alpha += (1. - T_i) * (1. - alpha);
}

float terrain_map(
	_in(vec3) pos
){
	float h = fbm(pos * 2., 2.);
	float hs = smoothstep(.35, 1., h);

	return hs;
}

vec3 render_planet(
	_in(ray_t) eye
){
	const float max_height = .25;

	sphere_t atmosphere = planet;
	atmosphere.radius += max_height;

	hit_t hit = no_hit;
	intersect_sphere(eye, atmosphere, hit);
	if (hit.material_id < 0) {
		return background(eye);
	}

#if 0 // test with checkboard pattern
	vec3 d = mul(rotate_around_x(u_time * 16.), hit.normal);
#ifdef HLSL
#define atan(y, x) atan2(x, y)
#endif
	float u = .5 + atan(d.z, d.x) / (2. * PI);
	float v = .5 - asin(d.y) / (1. * PI);
	float n = checkboard_pattern(vec2(u, v), 20.);
	return vec3(n, n, n);
#endif

	float t_min = .01;
	float t_max = max_height * 3;
	float t_step = t_max / 20.;

	float T = 1.;
	vec3 C = vec3(0, 0, 0);
	float alpha = 0.;

	for (float t = t_min; t < t_max; t += t_step) {
		vec3 p = hit.origin + t * eye.direction;

		vec3 n = p;// - planet.origin;
		normalize(n);
		n = mul(rotate_around_x(u_time * 16.), n);

		float hs = terrain_map(n);
		float h = planet.radius + hs * max_height;

		clouds_map(n, T, C, alpha, t_step);

		if (dot(p, p) < (h*h)) {
			vec3 terr = mix(
				vec3(.1, .1, .9),
				vec3(1, 1, 1),
				hs);
			return mix(terr, C, alpha);
		}

		//t_step += .001 * t;
	}

	return mix(background(eye), C, alpha);
}

vec3 render(
	_in(ray_t) eye
){
	return render_planet(eye);
}

#include "main.h"