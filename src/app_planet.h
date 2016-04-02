#include "def.h"
#include "util.h"
#include "intersect.h"

#include "noise_iq.h"
#define noise(x) noise_iq(x)
#include "fbm.h"

_constant(sphere_t) planet = _begin(sphere_t)
	vec3(0, 0, 0), 1., 0
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
	_in(float) t_step,
	_in(float) h
){
	//TODO: add offset and rot matrix to distrib
	// the noise better - right now seems to
	// follow the terrain
	float dens = fbm(pos * 2.234343, 1.76001263);

	const float coverage = .575675;
	dens *= step(coverage, dens);
	dens *= band(.2, .4, .6, exp(h) / 4.);

	const float absorbtion = 33.93434;
	float T_i = exp(-absorbtion * dens * t_step);
	T *= T_i;
	C += T  * (exp(h) / .055)
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

vec3 terrain_normal(
	_in(vec3) p
){
	float dt = 0.001;
	vec3 x = vec3(dt, 0, 0);
	vec3 y = vec3(0, dt, 0);
	vec3 z = vec3(0, 0, dt);
#define sdf terrain_map
	return normalize(vec3(
		sdf(p + x) - sdf(p - x),
		sdf(p + y) - sdf(p - y),
		sdf(p + z) - sdf(p - z)
	));
}

vec3 illuminate(
	_in(vec3) V,
	_in(vec3) L,
	_in(hit_t) hit,
	_in(vec3) color
){
	vec3 diffuse = max(0., dot(L, hit.normal))
		* color;

	vec3 H = normalize(L + V);
	vec3 specular = pow(
		max(0., dot(H, hit.normal)),
		50.) * vec3(1, 1, 1);

	return diffuse + specular;
}

vec3 render_planet(
	_in(ray_t) eye
){
	const float max_height = .35;
	const mat3 rot = rotate_around_x(u_time * 16.);

	sphere_t atmosphere = planet;
	atmosphere.radius += max_height;

	hit_t hit = no_hit;
	intersect_sphere(eye, atmosphere, hit);
	if (hit.material_id < 0) {
		return background(eye);
	}

#if 0 // test with checkboard pattern
	vec3 d = mul(rot, hit.normal);
#ifdef HLSL
#define atan(y, x) atan2(x, y)
#endif
	float u = .5 + atan(d.z, d.x) / (2. * PI);
	float v = .5 - asin(d.y) / (1. * PI);
	float n = checkboard_pattern(vec2(u, v), 20.);
	return vec3(n, n, n);
#endif

	const float t_min = .01;
	const float t_max = max_height * 3; //TODO: optimal value
	const float t_step = t_max / 120.; //TODO: optimal num of steps

	//TODO: better names (way to handle this)
	float T = 1.;
	vec3 C = vec3(0, 0, 0);
	float alpha = 0.;

	for (float t = t_min; t < t_max; t += t_step) {
		vec3 p = hit.origin + t * eye.direction;
		vec3 n = mul(rot, p);

		float hs = terrain_map(n);
		float h = planet.radius + hs * max_height;

		float p_len = length(n); //TODO: possible to get rid of?
		clouds_map(n, T, C, alpha, t_step,
			(p_len - planet.radius) / max_height);

		if (p_len < h) {
			//TODO: find more accurate intersection
			// A.linear interpolate from prev data
			// B. bsearch like https://www.shadertoy.com/view/4slGD4
			
			//TODO: fix normals
			//vec3 tn = terrain_normal(p);
			//if (dot(tn, tn) < .01*.01) {
			//	tn = n;
			//}
			
			hit_t impact = _begin(hit_t)
				t, // ray length at impact
				0, // material id
				n, // normal
				p // point of impact				
			_end;
#if 0			
			vec3 terr = illuminate(
				eye.direction,
				mul(rotate_around_y(u_time * 16.),
				vec3 (1, 0, 0)),
				impact,
				vec3 (1, 1, 1));
#endif

			vec3 snow = mix(
				vec3(.5, .5, .5),
				vec3(1, 1, 1),
				smoothstep(.75, 1., hs));
			vec3 rock = mix(
				vec3(1, 0, 0),
				snow,
				smoothstep(.4, .75, hs));			
			vec3 grass = mix(
				vec3(0, 1, 0),
				rock,
				smoothstep(.1, .4, hs));				
			vec3 c = mix(
				vec3(.1, .1, .9),
				grass,
				smoothstep (0, .1, hs));
				
			return mix(c, C, alpha);
		}

		//t_step += .001 * t; //TODO: research adaptive step/error
	}

	return mix(background(eye), C, alpha);
}

vec3 render(
	_in(ray_t) eye
){
	return render_planet(eye);
}

#include "main.h"