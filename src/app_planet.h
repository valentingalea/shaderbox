#include "def.h"
#include "util.h"
#include "intersect.h"

#include "noise_iq.h"
#define noise(x) noise_iq(x)
#include "fbm.h"

_constant(sphere_t) planet = _begin(sphere_t)
	vec3(0, 0, 0), 1., 0
_end;
_constant(float) max_height = .35;

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

struct cloud_drop_t {
	vec3 pos;
	float h;
	float T;
	vec3 C;
	float alpha;
};

cloud_drop_t start_cloud(
	_in(vec3) p
){
	return _begin(cloud_drop_t)
		p, 0., 1., vec3 (0., 0., 0.), 0.
	_end;
}

_mutable(cloud_drop_t) cloud;

void clouds_map(
	_inout(cloud_drop_t) cloud,
	_in(float) t_step
){
	float H = exp(cloud.h);
	float dens = fbm(cloud.pos * 1.2343 
		+ vec3(1.35, 3.35, 2.67), 2.9760);

	const float coverage = .575675; // higher=less clouds
	const float fuzziness = .0335; // higher=fuzzy, lower=blockier
	//dens *= step(coverage, dens);
	dens *= smoothstep(coverage, coverage + fuzziness, dens);

	// these 2 are identical ways to "band" vertically
	//TODO: understand why the exp thing really works
	//dens *= band(.2, .4, .6, exp(h) / 4.);
	dens *= 1. - smoothstep(.4, .6, H / 4.);

	const float absorbtion = 33.93434;
	float T_i = exp(-absorbtion * dens * t_step);
	cloud.T *= T_i;
	cloud.C +=
		cloud.T * (H / .055)
		* dens * t_step;
	cloud.alpha += (1. - T_i) * (1. - cloud.alpha);
}

float terrain_map(
	_in(vec3) pos
){
	float h = fbm(pos * 2., 2.);
	float hs = smoothstep(.35, 1., h);

	return hs;
}

float sdf_map(
	_in(vec3) pos
){
	float n = terrain_map(pos);
	return length(pos) - planet.radius - n * max_height;
}

vec3 terrain_normal(
	_in(vec3) p
){
#define F sdf_map
	vec3 dt = vec3(0.001, 0, 0);

	return normalize(vec3(
		F(p + dt.xzz) - F(p - dt.xzz),
		F(p + dt.zxz) - F(p - dt.zxz),
		F(p + dt.zzx) - F(p - dt.zzx)
	));
#undef F
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
	mat3 rot = rotate_around_x(u_time * 16.);
	mat3 rot2 = rotate_around_x(u_time * -8.);

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
	const float t_max = max_height * 4.; //TODO: optimal value
	const float t_step = t_max / 120.; //TODO: optimal num of steps
	vec3 prev_p = hit.origin;

	cloud = start_cloud (hit.origin);

	for (float t = t_min; t < t_max; t += t_step) {
		vec3 o = hit.origin + t * eye.direction;
		vec3 p = mul(rot, o - planet.origin);

		float hs = terrain_map(p);
		float h = planet.radius + hs * max_height;

		float p_len = length(p); //TODO: possible to get rid of?
#if 1
		cloud.pos = mul(rot2, o);
		cloud.h = (p_len - planet.radius) / max_height;
		clouds_map(cloud, t_step);
#endif

		if (p_len < h) {
			//TODO: find more accurate intersection
			// bsearch like https://www.shadertoy.com/view/4slGD4
			vec3 H = prev_p + (prev_p - p) * .5;
			float hs = terrain_map(H);
			
			vec3 n = terrain_normal(H);
			//return n;
			
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
				smoothstep (.0, .1, hs));
#if 0
			hit_t impact = _begin(hit_t)
				t, // ray length at impact
				0, // material id
				n, // normal
				p // point of impact				
			_end;
			
			c = illuminate(
				eye.direction,
				vec3(1, 0, 0),
				impact,
				c);
#endif
			return mix(c, cloud.C, cloud.alpha);
		}

		prev_p = p;
		//t_step += .001 * t; //TODO: research adaptive step/error
	}

	return mix(background(eye), cloud.C, cloud.alpha);
}

vec3 render(
	_in(ray_t) eye
){
	return render_planet(eye);
}

#include "main.h"