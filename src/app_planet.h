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

vec3 render(
	_in(ray_t) eye
){
	//return background(eye);

	hit_t hit = no_hit;
	intersect_sphere(eye, planet, hit);
	if (hit.material_id < 0) {
		return background (eye);
	}
	
	vec3 d = mul(rotate_around_x(u_time * 16.), hit.normal);
#ifdef HLSL
#define atan(y, x) atan2(x, y)
#endif
	float u = .5 + atan(d.z, d.x) / (2. * PI);
	float v = .5 - asin(d.y) / (1. * PI);
	float n = //checkboard_pattern(vec2(u, v), 50.);
	fbm (vec3 (u, v, 0) * 4., 3.);
	
	float s = smoothstep (.45, .5, n);
	vec3 color = mix (
		vec3(.1, .1, .9),
		vec3(n),
		s);
	
	return color;
}

#include "main.h"