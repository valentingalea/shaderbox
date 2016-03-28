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
	eye = vec3 (0, 0, -2.5);
	look_at = vec3 (0, 0, 2);
}

float density_func(
	_in(vec3) pos,
//	_in(vec3) offset,
	_in(float) coverage,
	_in(float) fuziness
){
	vec3 p = pos;// * .0212242 + offset;
	float dens = fbm(p * 6., 1.);
	
	//dens *= step(coverage, dens);
	//dens -= coverage;
	dens *= smoothstep (coverage, coverage + fuziness, dens);

	return dens;// clamp(dens, 0., 1.);	
}

vec4 render_clouds(
	_in(ray_t) eye,
	_in(float) coverage,
	_in(float) thickness,
	_in(float) absorbtion,
	_in(float) fuzziness
){
	hit_t hit = no_hit;
	sphere_t atmosphere = planet;
	atmosphere.radius += .25;
	
	intersect_sphere(eye, atmosphere, hit);
	if (hit.material_id < 0) {
		return vec4 (0, 0, 0, 0);
	}

	const int steps = 5;
	float march_step = thickness / float(steps);
	vec3 dir_step = eye.direction * march_step;
	vec3 pos = hit.origin;

	float T = 1.;
	vec3 C = vec3(0, 0, 0);
	float alpha = 0.;

	for (int i = 0; i < steps; i++) {
		float h = float(i) / float(steps);
		float dens = density_func(
			pos, coverage, fuzziness);

		float T_i = exp(-absorbtion * dens * march_step);
		T *= T_i;
		//if (T < .01) break;

		C += T * 
		//	(exp(h) / 1.75) * // fake light
			dens * march_step;
		alpha += (1. - T_i) * (1. - alpha);

		pos += dir_step;
		//if (length(pos) > 1e3) break;
	}

	return vec4(C, alpha);
}

vec3 render_planet(
	_in(ray_t) eye
){
	const float max_height = .2;
	
	hit_t hit = no_hit;
	sphere_t atmosphere = planet;
	atmosphere.radius += max_height;
	intersect_sphere(eye, atmosphere, hit);
	if (hit.material_id < 0) {
		return background (eye);
	}
/*	
	vec3 d = mul(rotate_around_x(u_time * 16.), hit.normal);
#ifdef HLSL
#define atan(y, x) atan2(x, y)
#endif
#if 0
	float u = .5 + atan(d.z, d.x) / (2. * PI);
	float v = .5 - asin(d.y) / (1. * PI);
	float n = checkboard_pattern(vec2(u, v), 20.);
	vec3 color = vec3 (n, n, n);
#else
	float n = fbm (d * 4., 2.);
	
	float s = smoothstep (.45, .5, n);
	vec3 color = mix (
		vec3(.1, .1, .9),
		vec3(n, n, n),
		s);
#endif
*/
	float t_min = .01;
	float t_max = max_height * 3;
	float t_step = t_max / 20.;
	
	float T = 1.;
	vec3 C = vec3(0, 0, 0);
	float alpha = 0.;
	
	for (float t = t_min; t < t_max; t += t_step) {
		vec3 p = hit.origin + t * eye.direction;
	
		vec3 n = p;// - planet.origin;
		normalize (n);
		n = mul(rotate_around_x(u_time * 4.), n);
		
		float h = fbm (n * 2., 2.);
		float hs = smoothstep (.35, 1., h);
		float hhs = hs * max_height;
		h = planet.radius + hhs;
		
		float dens = density_func(n, .75634, .035);
		float T_i = exp(-5.03 * dens * t_step);
		T *= T_i;
		C += T * (exp(hs) / 1.75) * // fake light
			dens * t_step;
		alpha += (1. - T_i) * (1. - alpha);
		
		if (dot (p, p) < (h*h)) {
			vec3 a =  mix (
				vec3(.1, .1, .9),
				vec3(1, 0, 0),
				hs);
			return mix (a, C, alpha);
		}
		//t_step += .001 * t;
	}
	return mix (background (eye), C, alpha);
}

float terrain (vec2 p)
{
	vec3 d = vec3 (p.x, 0, p.y);
	float n = fbm (d * 2. + u_time, 2.);
	float s = smoothstep (.45, .5, n);
	return s;
}

vec3 render(
	_in(ray_t) eye
){/*
	float t_min = .1;
	float t_max = 5.;
	const float t_step = .125;
	for (float t = t_min;
	t < t_max; t += t_step) {
		vec3 p = eye.origin + t * eye.direction;
		float h = terrain (p.xz);
		if (p.y < h) {
			return vec3 (h);
		}
	}
	return vec3 (0, .56, 0);
*/
	vec3 pln;
	return render_planet(eye);

	vec4 cld = render_clouds(
		eye,
		1. - .55,
		.25,
		11.0345,
		.035
	);
	
	return mix(pln, cld.rgb, cld.a);
}

#include "main.h"