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
	look_at = vec3 (0, 0, -1);
}

float density_func(
	_in(vec3) pos,
//	_in(vec3) offset,
	_in(float) coverage,
	_in(float) fuziness
){
	vec3 p = pos;// * .0212242 + offset;
	float dens = fbm(p * 4., 3);
	
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
	hit_t hit = no_hit;
	intersect_sphere(eye, planet, hit);
	if (hit.material_id < 0) {
		return background (eye);
	}
	
	vec3 d = mul(rotate_around_x(u_time * 16.), hit.normal);
#ifdef HLSL
#define atan(y, x) atan2(x, y)
#endif
#if 1
	float u = .5 + atan(d.z, d.x) / (2. * PI);
	float v = .5 - asin(d.y) / (1. * PI);
	float n = checkboard_pattern(vec2(u, v), 20.);
	vec3 color = vec3 (n, n, n);
#else
	float n = fbm (d * 4., 3.);
	
	float s = smoothstep (.45, .5, n);
	vec3 color = mix (
		vec3(.1, .1, .9),
		vec3(n, n, n),
		s);
#endif
	
	return color;
}

vec3 render(
	_in(ray_t) eye
){
	vec3 pln = render_planet(eye);

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