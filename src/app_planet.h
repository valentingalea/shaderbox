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
	float t;
	vec3 origin;
	vec3 pos;
	float h;
	float T;
	vec3 C;
	float alpha;
};

cloud_drop_t start_cloud(
	_in(vec3) origin
){
	cloud_drop_t c = _begin(cloud_drop_t)
		0., origin, origin, 0., 1., vec3 (0., 0., 0.), 0.
	_end;
	return c;
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
	dens *= 1. - smoothstep(.4, .5, H / 4.);

	const float absorbtion = 33.93434;
	float T_i = exp(-absorbtion * dens * t_step);
	cloud.T *= T_i;
	cloud.C +=
		cloud.T * (H / .055)
		* dens * t_step;
	cloud.alpha += (1. - T_i) * (1. - cloud.alpha);
}

void clouds_march(
	_in(ray_t) eye,
	_inout(cloud_drop_t) cloud,
	_in(float) exit_dist
){
	const int steps = 75;
	float t_step = (max_height * 4.) / float(steps);
	mat3 rot2 = rotate_around_x(u_time * -8.);

	for (int i = 0; i < steps; i++) {
		vec3 o = cloud.origin + cloud.t * eye.direction;
		cloud.pos = mul(rot2, o - planet.origin);

		cloud.h = (length(cloud.pos) - planet.radius) / max_height;
		cloud.t += t_step;

		clouds_map(cloud, t_step);
		if (cloud.t > exit_dist) return;
	}
}

float terrain_map(
	_in(vec3) pos
){
	float h = fbm(pos * 2., 2.);
	float hs = smoothstep(.35, 1., h);

	return hs;
}

vec2 sdf_map(
	_in(vec3) pos
){
	float n = terrain_map(pos);
	return vec2(length(pos) - planet.radius - n * max_height, n);
}

vec3 terrain_normal(
	_in(vec3) p
){
#define F(t) sdf_map(t).x
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

	cloud = start_cloud (hit.origin);

	float t = 0.;
	for (int i = 0; i < 30; i++) {
		vec3 o = hit.origin + t * eye.direction;
		vec3 p = mul(rot, o - planet.origin);

		vec2 d = sdf_map(p);

		clouds_march(eye, cloud, max_height);

		if (d.x < .01) {
			float hs = d.y;
			vec3 n = terrain_normal(p);
			//return n;

			clouds_march(eye, cloud, t);
			
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

		t += d.x;
	}

	clouds_march(eye, cloud, max_dist);
	return mix(background(eye), cloud.C, cloud.alpha);
}

vec3 render(
	_in(ray_t) eye,
	_in(vec3) point_cam
){
	return render_planet(eye);
}

#define FOV tan(radians(30.))
#include "main.h"