#include "def.h"
#include "util.h"
#include "noise.h"

#define u_time iGlobalTime

vec3 background(_in(ray_t) ray)
{
	return vec3(.1, .1, .7);
}

void setup_scene()
{
}

void setup_camera(_inout(vec3) eye, _inout(vec3) look_at)
{
	mat3 rot = rotate_around_y (u_time * 10.);
	eye = rot * vec3(0, 10, 50);
	look_at = vec3(0, 0, 0);
}

vec3 illuminate(_in(hit_t) hit)
{
	vec3 L = normalize (vec3 (1, 1, 0));
	return max (0, dot (hit.normal, L))
	* vec3 (0, 1, 0);
}

float terrain_func(_in(vec2) t)
{
	return
	//NOISE (t/3.)
	sin (t.x) * sin (t.y)
	;
}

vec3 terrain_normal(_in(vec2) p)
{
#define EPS 0.05
	vec2 dx = vec2 (EPS, 0);
	vec2 dz = vec2 (0, EPS);
	
	vec3 n = vec3 (
		terrain_func (p - dx) - terrain_func(p + dx),
		2. * EPS,
		terrain_func (p - dz) - terrain_func(p + dz)
	);
	
	return normalize  (n);
}

vec3 render(_in(ray_t) ray)
{
	const float min_t = .01;
	const float max_t = 50.;
	float step = .5;

	for (float t = min_t; t < max_t; t += step) {
		vec3 p = ray.origin + ray.direction * t;
		float h = terrain_func (p.xz);
	
		if (p.y < h) {
			hit_t h = hit_t _begin
				t, // ray length at impact
				0, // material id
				float(t) / float(max_t), // material custom param
				terrain_normal(p.xz),
				p // point of impact				
			_end;
			
			return illuminate  (h);
		}
		
		// decay accuracy inversly prop with distance
		//step += h;
	}

	return background(ray);
}

#include "main.h"