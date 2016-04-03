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
	float dens = fbm(pos * 1.2343 
		+ vec3(1.35, 3.35, 2.67), 2.9760);

	const float coverage = .575675;
	dens *= step(coverage, dens);
	//dens *= smoothstep(coverage, coverage + 0.0343, dens);
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
#define F terrain_map
#if 1
	vec3 dt = vec3(0.001, 0, 0);

	return normalize(vec3(
		F(p + dt.xzz) - F(p - dt.xzz),
		F(p + dt.zxz) - F(p - dt.zxz),
		F(p + dt.zzx) - F(p - dt.zzx)
	));
#else
	const float dt = 0.001;

	vec3 n = normalize(p);
	vec3 tangent, binormal;
	fast_orthonormal_basis(n, tangent, binormal);

	return normalize(vec3(
		F(p + binormal*dt) - F(p - binormal*dt),
		F(p + n*dt) - F(p - n*dt),
		F(p + tangent*dt) - F(p - tangent*dt)
	));
#endif
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
	const float max_height = .35;
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

	//TODO: better names (way to handle this)
	float T = 1.;
	vec3 C = vec3(0, 0, 0);
	float alpha = 0.;

	for (float t = t_min; t < t_max; t += t_step) {
		vec3 o = hit.origin + t * eye.direction;
		vec3 p = mul(rot, o - planet.origin);

		float hs = terrain_map(p);
		float h = planet.radius + hs * max_height;

		float p_len = length(p); //TODO: possible to get rid of?
#if 1
		clouds_map(
			mul(rot2, o),
			T, C, alpha, t_step,
			(p_len - planet.radius) / max_height);
#endif

		if (p_len < h) {
			//TODO: find more accurate intersection
			// A.linear interpolate from prev data
			// B. bsearch like https://www.shadertoy.com/view/4slGD4
			vec3 H = prev_p + (prev_p - p) * .5;
			float hs = terrain_map(H);
			
			//TODO: fix normals
			//vec3 n = terrain_normal(H);
			//if (dot(tn, tn) < .01*.01) tn = n;
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
			return mix(c, C, alpha);
		}

		prev_p = p;
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