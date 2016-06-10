#define APP_CLOUDS
#include "def.h"
#include "util.h"
#include "intersect.h"

#define hg_g (.2)
#include "volumetric.h"

//#define SKY_SPHERE
#define USE_NOISE_TEX

// ----------------------------------------------------------------------------
// Scene
// ----------------------------------------------------------------------------
#ifdef SKY_SPHERE
_constant(sphere_t) atmosphere = _begin(sphere_t)
	vec3(0, -499, 0), 500., 0
_end;
#define cld_noise_factor (1. / atmosphere.radius)
#else
#define cld_noise_factor .001
#endif

void setup_camera(
	_inout(vec3) eye,
	_inout(vec3) look_at
){
	eye = vec3(0, -.5, 0);
	float angle = u_mouse.x * .5;
	look_at = mul(rotate_around_y(angle), vec3(0, 0, -1));
}

void setup_scene()
{
}

vec3 render_sky_color(
	_in(vec3) eye_dir
){
	float sun_amount = max(dot(eye_dir, sun_dir), 0.);

	vec3 sky = mix(vec3(.0, .1, .4), vec3(.3, .6, .8), 1.0 - eye_dir.y);
	sky += sun_color * min(pow(sun_amount, 1500.0) * 5.0, 1.0);
	sky += sun_color * min(pow(sun_amount, 10.0) * .6, 1.0);

	return abs(sky);
}

// ----------------------------------------------------------------------------
// Density
// ----------------------------------------------------------------------------
#ifdef USE_NOISE_TEX
Texture3D u_tex_noise : register(t1);
Texture3D u_tex_noise_2 : register(t2);
SamplerState u_sampler0 : register(s0);
#else
#include "noise_iq.h"
#include "fbm.h"
DECL_FBM_FUNC(fbm, 4, noise_iq(p))
#endif

float density_func(
	_in(vec3) pos_in,
	_in(float) height
){
	vec3 pos = pos_in * cld_noise_factor - wind_dir * u_time;

	float shape =
#ifdef USE_NOISE_TEX
	u_tex_noise.SampleLevel(u_sampler0, pos * .7015460, 0).r;
#else
	fbm(pos * 2.03, 2.64, .5, .5);
#endif

#ifdef USE_NOISE_TEX
	float w =
		//fbm_worley_tile(pos, 7., 1., .5);
		u_tex_noise_2.SampleLevel(u_sampler0, pos * 1.73547, 0).r;
	float ww = mix(w, 1. - w, height);// exp(height) / 3.23);
	shape = remap(shape, ww * .7, 1., 0., 1.);
#endif

	return shape * smoothstep(cld_coverage, cld_coverage + .0135, shape);
	//return smoothstep(cld_coverage, 1., shape);
}

// ----------------------------------------------------------------------------
// Volumetrics
// ----------------------------------------------------------------------------
float illuminate_volume(
	_in(vec3) origin,
	_in(float) height,
	_in(vec3) V,
	_in(vec3) L
){
#if 0
	float luminance = exp(height) / 2.;
#else
	const float dt = cld_thick / float(cld_march_steps);
	volume_sampler_t vol = construct_volume(origin);
	vol.pos += L * dt; // don't sample just where the main raymarcher is

#ifdef HLSL
	[fastopt] [loop]
#endif
	for (int i = 0; i < illum_march_steps; i++) {
		vol.height = float(i) / float(illum_march_steps);
		float density = density_func(vol.pos, vol.height);

		vol.transmittance *= exp(- density * sigma_scattering * dt);
		vol.pos += L * dt;
	}

	float luminance = vol.transmittance;
#endif

#if 0
	return luminance;
#else
	return luminance * sun_power * henyey_greenstein_phase_func(clamp(dot(L, V), 0., 1.));
#endif
}

void integrate_volume(
	_inout(volume_sampler_t) vol,
	_in(vec3) V,
	_in(vec3) L,
	_in(float) density,
	_in(float) dt
){
	if (density < .005) return;

	// change in transmittance (follows Beer-Lambert law)
	float T_i = exp(-density * sigma_scattering * dt);
	// Update accumulated transmittance
	vol.transmittance *= T_i;

	// integrate output radiance (here essentially color)
	vol.radiance += 
		(density * sigma_scattering) * 
		illuminate_volume(vol.pos, vol.height, V, L) *
		vol.transmittance * 
		dt;

	// accumulate opacity
	vol.alpha += (1. - T_i) * (1. - vol.alpha);
}

// ----------------------------------------------------------------------------
// Raymarching
// ----------------------------------------------------------------------------
vec4 render_clouds(
	_in(ray_t) eye
){
	const float march_step = cld_thick / float(cld_march_steps);

#ifdef SKY_SPHERE
	hit_t hit = no_hit;
	intersect_sphere_from_inside(eye, atmosphere, hit);

	vec3 projection = eye.direction;
	vec3 origin = hit.origin;
#else
	vec3 projection = eye.direction / eye.direction.y;
	vec3 origin = eye.origin + projection * 150.;
#endif

	vec3 iter = projection * march_step;
	volume_sampler_t cloud = construct_volume(origin);

#ifdef HLSL
	[fastopt] [loop]
#endif
	for (int i = 0; i < cld_march_steps; i++) {
		cloud.height = float(i) / float(cld_march_steps);
		
		float density = density_func(cloud.pos, cloud.height);

		integrate_volume(
			cloud,
			eye.direction,
			sun_dir,
			density,
			march_step);

		cloud.pos += iter;

		if (cloud.alpha > .999) break;
	}

	float cutoff = dot(eye.direction, vec3(0, 1, 0));
	return vec4(cloud.radiance, cloud.alpha * smoothstep(.0, .2, cutoff));
}

vec3 render(
	_in(ray_t) eye_ray,
	_in(vec3) point_cam
){
	vec3 sky = render_sky_color(eye_ray.direction);
#ifdef HLSL
	[flatten]
#endif
	if (dot(eye_ray.direction, vec3(0, 1, 0)) < 0.05) return sky;

	vec4 cld = render_clouds(eye_ray);
	vec3 col = mix(sky, cld.rgb, cld.a);

	return abs(col);
}

#define FOV tan(radians(30.))
#include "main.h"