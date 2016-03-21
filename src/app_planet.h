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
	_in(ray_t) ray
){
	return vec3(0, 1, 0);
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
	hit_t hit = no_hit;
	intersect_sphere(eye, planet, hit);
	if (hit.material_id < 0) {
		return background (eye);
	}
	
	float n = fbm (hit.origin, 4);
	vec3 color = vec3(n);
	
	return color;
}

#include "main.h"