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
#include "volumetric.h"

#ifdef NOISE_VALUE
#include "noise_iq.h"
#define noise(x) noise_iq(x)
#endif
#ifdef NOISE_WORLEY
#include "noise_worley.h"
#define noise(x) (1. - noise_w(x).r)
//#define noise(x) abs( noise_iq(x / 8.) - (1. - (noise_w(x * 2.).r)))
#endif

#include "fbm.h"
DECL_FBM_FUNC(fbm_clouds, 4, .5)

#ifdef HLSLTOY
Texture3D u_tex_noise : register(t1);
SamplerState u_sampler0 : register(s0);
#endif

float noise_func(_in(vec3) x)
{
#if 0
	return u_tex_noise.Sample(u_sampler0, x);
#else
	return fbm_clouds(x, FBM_FREQ);
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

	return dens;
}

float illuminate_volume(
	_inout(volume_sampler_t) cloud,
	_in(vec3) V,
	_in(vec3) L
){
	//return exp(cloud.height) / 1.75;
	
	const int steps = 5;
	const float t_step = 1.;
	float t = 0.;
	float R = 0.;
	
	volume_sampler_t vol = begin_volume(
		cloud.pos, ABSORPTION);
	
	float mu = dot (V, L);
	float phase = henyey_greenstein_phase_func (mu);

	for (int i = 0; i < steps; i++) {
		vec3 p = vol.origin + t * L;
		float dens = density_func(p, WIND_DIR, COVERAGE, FUZZINESS);

		// change in transmittance (follows Beer-Lambert law)
		float T_i = exp(-vol.coeff_absorb * dens * t_step);
		// Update accumulated transmittance
		vol.T *= T_i;
		// integrate output radiance
		R += vol.T * vol.coeff_absorb * phase;
	
		t += t_step;
	}
	
	return R;
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

	volume_sampler_t cloud = begin_volume(pos, ABSORPTION);

	for (int i = 0; i < steps; i++) {
		float dens = density_func(cloud.pos, wind_dir, coverage, fuzziness);

		cloud.height = float(i) / float(steps);
		integrate_volume(cloud,
			eye.direction, sun_dir,
			dens, march_step);

		cloud.pos += dir_step;
	}

	return vec4(cloud.C, cloud.alpha);
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