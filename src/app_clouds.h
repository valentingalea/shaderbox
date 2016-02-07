/**** TWEAK *****************************************************************/
#define COVERAGE		.50
#define THICKNESS		15.
#define ABSORPTION		1.030725
#define WIND_DIR		vec3(0, 0, -u_time * .2)

#define FBM_FREQ		2.76434
#define NOISE_VALUE
//#define NOISE_WORLEY
//#define NOISE_PERLIN

//#define SIMULATE_LIGHT
#define FAKE_LIGHT
#define SUN_DIR			normalize(vec3(0, abs(sin(u_time * .3)), -1))

#define STEPS			25
/******************************************************************************/

#include "def.h"
#include "util.h"
#include "intersect.h"

#include "noise_iq.h"
//#include "noise_worley.h"
//#include "../lib/ashima-noise/src/classicnoise3D.glsl"

#ifdef NOISE_VALUE
#define noise(x) noise_iq(x)
#endif
#ifdef NOISE_WORLEY
#define noise(x) (1. - noise_w(x).r)
//#define noise(x) abs( noise_iq(x / 8.) - (1. - (noise_w(x * 2.).r)))
#endif
#ifdef NOISE_PERLIN
#define noise(x) abs(cnoise(x))
#endif
#include "fbm.h"

#ifdef HLSLTOY
Texture3D u_tex_noise : register(t1);
SamplerState u_sampler0 : register(s0);
#endif

float noise_func(_in(vec3) x)
{
#if 0
	return u_tex_noise.Sample(u_sampler0, x);
#else
	return fbm(x, FBM_FREQ);
#endif
}

vec3 render_sky_color(
	_in(vec3) eye_dir,
	_in(vec3) sun_dir
){
	_constant(vec3) sun_color = vec3(1., .7, .55);

	float sun_amount = max(dot(eye_dir, sun_dir), 0.);

	vec3  sky = mix(vec3(.0, .1, .4), vec3(.3, .6, .8), 1.0 - eye_dir.y);
	sky = sky + sun_color * min(pow(sun_amount, 1500.0) * 5.0, 1.0);
	sky = sky + sun_color * min(pow(sun_amount, 10.0) * .6, 1.0);

	return sky;
}

float density_func(
	_in(vec3) pos,
	_in(vec3) offset,
	_in(float) coverage
){
	vec3 p = pos * .0212242 + offset;
	float dens = noise_func(p);
	
	//dens = band (.1, .3, .6, dens);
	//dens *= step(cov, dens);
	//dens -= cov;
	dens *= smoothstep (coverage, coverage + .05, dens);

	return clamp(dens, 0., 1.);	
}

float gather_light(
	_in(vec3) origin,
	_in(vec3) sun_dir,
	_in(vec3) wind_dir,
	_in(float) coverage,
	_in(float) absorbtion
){
	const int steps = 8;
	const float march_step = 1.;

	vec3 pos = origin;
	vec3 dir_step = sun_dir * march_step;
	float T = 1.;

	for (int i = 0; i < steps; i++) {
		float dens = density_func(pos, wind_dir, coverage);

		float T_i = exp(-absorbtion * dens * march_step);
		T *= T_i;

		pos += dir_step;
	}

	return T;
}

vec4 render_clouds(
	_in(ray_t) eye,
	_in(vec3) sun_dir,
	_in(vec3) wind_dir,
	_in(float) coverage,
	_in(float) thickness,
	_in(float) absorbtion
){
#if 0 // atmosphere 'sphere' intersect
	_constant(sphere_t) atmosphere = _begin(sphere_t)
		vec3(0, -450, 0), 500., 0
	_end;
	_constant(sphere_t) atmosphere_2 = _begin(sphere_t)
		atmosphere.origin, atmosphere.radius + float(50.), 0
	_end;

	hit_t hit = no_hit;
	intersect_sphere(eye, atmosphere, hit);
	hit_t hit_2 = no_hit;
	intersect_sphere(eye, atmosphere_2, hit_2);

	const float thickness = length(hit_2.origin - hit.origin);
	const float r = 1. - ((atmosphere_2.radius - atmosphere.radius) / thickness);
	const int steps = STEPS + int(32. * r);
	float march_step = thickness / float(steps);
	vec3 dir_step = eye.direction * march_step;
	vec3 pos = hit.origin;
#else // plane projection
	const int steps = STEPS;
	float march_step = thickness / float(steps);
	vec3 dir_step = eye.direction / eye.direction.y * march_step;
	vec3 pos = eye.origin + eye.direction * 100.;
#endif

	float T = 1.;
	vec3 C = vec3(0, 0, 0);
	float alpha = 0.;

	for (int i = 0; i < steps; i++) {
		float h = float(i) / float(steps);
		float dens = density_func (pos, wind_dir, coverage);

		float T_i = exp(-absorbtion * dens * march_step);
		T *= T_i;
		//if (T < .01) break;

		C += T * 
#ifdef SIMULATE_LIGHT
			gather_light(pos) *
#endif
#ifdef FAKE_LIGHT
			(exp(h) / 1.75) *
#endif
			dens * march_step;
		alpha += (1. - T_i) * (1. - alpha);

		pos += dir_step;
		//if (length(pos) > 1e3) break;
	}

	return vec4(C, alpha);
}

#ifdef UE4
vec3 ue4_render_clouds(
	_in(vec3) cam_dir,
	_in(float) time,
	_in(float) coverage,
	_in(float) thickness,
	_in(float) absorbtion,
	_in(vec3) sun_dir,
	_in(vec3) wind_dir
){
	ray_t eye_ray = _begin(ray_t)
		vec3(0, 0, 0),
		cam_dir
	_end;
	u_time = time;

	vec3 sky = render_sky_color(
		eye_ray.direction,
		sun_dir
		);
	vec4 cld = render_clouds(
		eye_ray,
		sun_dir,
		wind_dir,
		1. - coverage,
		thickness,
		absorbtion
	);
	return mix(sky, cld.rgb, cld.a);
}
#endif

void mainImage(
	_out(vec4) fragColor,
#ifdef SHADERTOY
	vec2 fragCoord
#else
	_in(vec2) fragCoord
#endif
){
#if 0 // DEBUG
	float n = saturate(get_noise(point_cam));
	fragColor = vec4(vec3(n, n, n), 1);
	return;
#endif

	vec3 col = vec3(0, 0, 0);

	vec2 aspect_ratio = vec2(u_res.x / u_res.y, 1);
	float fov = tan(radians(45.0));
	vec2 point_ndc = fragCoord.xy / u_res.xy;
#ifdef HLSLTOY
	point_ndc.y = 1. - point_ndc.y;
#endif
	vec3 point_cam = vec3((2.0 * point_ndc - 1.0) * aspect_ratio * fov, -1.0);

	vec3 eye = vec3(0, 1., 0);
	vec3 look_at = vec3(0, 1.6, -1);
	ray_t eye_ray = get_primary_ray(point_cam, eye, look_at);

	eye_ray.direction.yz = mul(rotate_2d(+u_mouse.y * .13), eye_ray.direction.yz);
	eye_ray.direction.xz = mul(rotate_2d(-u_mouse.x * .33), eye_ray.direction.xz);

	hit_t hit = no_hit;
	_constant(plane_t) ground = _begin(plane_t)
		vec3(0., -1., 0.), 0., 1
	_end;
	intersect_plane(eye_ray, ground, hit);

	if (hit.material_id == 1) {
		float cb = checkboard_pattern(hit.origin.xz, .5);
		col = mix(vec3(.6, .6, .6), vec3(.75, .75, .75), cb);
	} else {
#if 1
		vec3 sky = render_sky_color(
			eye_ray.direction,
			SUN_DIR
		);
		vec4 cld = render_clouds(
			eye_ray, 
			SUN_DIR,
			WIND_DIR,
			1. - COVERAGE,
			THICKNESS,
			ABSORPTION
		);
		col = mix(sky, cld.rgb, cld.a);
#else
		intersect_sphere(eye_ray, atmosphere, hit);
		vec3 d = hit.normal;
		float u = .5 + atan(d.z, d.x) / (2. * PI);
		float v = .5 - asin(d.y) / PI;
		float cb = checkboard_pattern(vec2(u, v), 50.);
		col = vec3(cb, cb, cb);
#endif
	}

	fragColor = vec4(col, 1);
}
