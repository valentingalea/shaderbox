/**** TWEAK *****************************************************************/
#define COVERAGE		.50
#define THICKNESS		15.
#define ABSORPTION		1.030725
#define FUZZINESS		0.035

#define FBM_FREQ		2.76434
#define NOISE_VALUE
//#define NOISE_WORLEY

#define WIND_DIR		vec3(0, 0, -u_time * .2)
#define SUN_DIR			normalize(vec3(0, abs(sin(u_time * .3)), -1))

#define STEPS			25
/******************************************************************************/

#include "def.h"
#include "util.h"
#include "intersect.h"

#include "noise_iq.h"
//#include "noise_worley.h"
#ifdef NOISE_VALUE
#define noise(x) noise_iq(x)
#endif
#ifdef NOISE_WORLEY
#define noise(x) (1. - noise_w(x).r)
//#define noise(x) abs( noise_iq(x / 8.) - (1. - (noise_w(x * 2.).r)))
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

	vec3 sky = mix(vec3(.0, .1, .4), vec3(.3, .6, .8), 1.0 - eye_dir.y);
	sky += sun_color * min(pow(sun_amount, 1500.0) * 5.0, 1.0);
	sky += sun_color * min(pow(sun_amount, 10.0) * .6, 1.0);

	return sky;
}

float density_func(
	_in(vec3) pos,
	_in(vec3) offset,
	_in(float) coverage,
	_in(float) fuziness
){
	vec3 p = pos * .0212242 + offset;
	float dens = noise_func(p);
	
	//dens *= step(coverage, dens);
	//dens -= coverage;
	dens *= smoothstep (coverage, coverage + fuziness, dens);

	return clamp(dens, 0., 1.);	
}

vec4 render_clouds(
	_in(ray_t) eye,
	_in(vec3) sun_dir,
	_in(vec3) wind_dir,
	_in(float) coverage,
	_in(float) thickness,
	_in(float) absorbtion,
	_in(float) fuzziness
){
#if 0 // atmosphere 'sphere' intersect
	_constant(sphere_t) atmosphere = _begin(sphere_t)
		vec3(0, -450, 0), 500., 0
	_end;

	hit_t hit = no_hit;
	intersect_sphere(eye, atmosphere, hit);

	const int steps = STEPS;
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
		float dens = density_func(pos, wind_dir, coverage, fuzziness);

		float T_i = exp(-absorbtion * dens * march_step);
		T *= T_i;
		//if (T < .01) break;

		C += T * 
			(exp(h) / 1.75) * // fake light
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
	_in(float) fuzziness,
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
		absorbtion,
		fuzziness
	);
	return mix(sky, cld.rgb, cld.a);
}
#endif

void setup_camera(
	_inout(vec3) eye,
	_inout(vec3) look_at
){
	eye = vec3(0, 1., 0);
	look_at = vec3(0, 1.6, -1);
}

void setup_scene()
{
	//eye_ray.direction.yz = mul(rotate_2d(+u_mouse.y * .13), eye_ray.direction.yz);
	//eye_ray.direction.xz = mul(rotate_2d(-u_mouse.x * .33), eye_ray.direction.xz);
}

vec3 render(
	_in(ray_t) eye_ray,
	_in(vec3) point_cam
){
	vec3 col = vec3(0, 0, 0);

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
			ABSORPTION,
			FUZZINESS
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

	return col;
}

#define FOV 1. // 45 degrees
#include "main.h"
