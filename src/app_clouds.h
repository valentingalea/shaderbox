#include "def.h"
#include "util.h"

#define hg_g (.2)
#include "volumetric.h"

// ----------------------------------------------------------------------------
// Scene
// ----------------------------------------------------------------------------
#define wind_dir	vec3(0, 0, -u_time * .2)
#define sun_dir		normalize(vec3(0, abs(sin(u_time * .3)), -1))
#define sun_color	vec3(1., .7, .55)

#define SPHERE
#ifdef SPHERE
_constant(sphere_t) atmosphere = _begin(sphere_t)
	vec3(0, -495, 0), 500., 0
_end;

void intersect_sphere(
	_in(ray_t) ray,
	_in(sphere_t) sphere,
	_inout(hit_t) hit
){
	vec3 rc = sphere.origin - ray.origin;
	float radius2 = sphere.radius * sphere.radius;
	float tca = dot(rc, ray.direction);
	float d2 = dot(rc, rc) - tca * tca;
	float thc = sqrt(radius2 - d2);
	float t0 = tca - thc;
	float t1 = tca + thc;

	vec3 impact = ray.origin + ray.direction * t0;
	hit.t = t0;
	hit.material_id = sphere.material;
	hit.origin = impact;
	hit.normal = (impact - sphere.origin) / sphere.radius;
}

#define cld_noise_factor (1. / atmosphere.radius)
#else
#define cld_noise_factor .001
#endif

void setup_camera(
	_inout(vec3) eye,
	_inout(vec3) look_at
){
	eye = vec3(0, -.5, 0);
	look_at = mul(rotate_around_y(u_mouse.x), vec3(0, 0, -1));
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
#define TEX
#ifdef TEX
Texture3D u_tex_noise : register(t1);
Texture3D u_tex_noise_2 : register(t2);
SamplerState u_sampler0 : register(s0);
#endif

#include "../lib/ashima-noise/src/common.glsl"
//#include "../lib/ashima-noise/src/classicnoise3d.glsl"
#include "../lib/ashima-noise/src/noise3d.glsl"
//#include "../lib/ashima-noise/src/cellular3d.glsl"
#include "noise_worley.h"

#include "fbm.h"
DECL_FBM_FUNC(fbm_simplex, 5, abs(snoise(p)))
DECL_FBM_FUNC_TILE(fbm_worley_tile, 4, (1. - (noise_w(p, L).r + .25)))


float density_func(
	_in(vec3) pos_in,
	_in(float) height
){
	vec3 pos = pos_in * cld_noise_factor;// -wind_dir;

	float base =
#ifdef TEX
	u_tex_noise.SampleLevel(u_sampler0, pos, 0).r;
#else
	fbm_simplex(pos * 2.03, 2.64, .5, .5);

	//float p = fbm_low_freq_perlin(pos * 4., 2., .5, .5);
	//float w = 1. - fbm_low_freq_worley(pos * 4., 4., .5, .5);
	//base = remap(p, -w, 1., 0., 1.);
#endif

#define cld_coverage (.735)
// my old method
	return smoothstep(cld_coverage, cld_coverage + .0135, base);
// book equiv method
	//return smoothstep(cld_coverage, 1., base);
// GPU Pro 7
	//float base_with_coverage = remap(base, cld_coverage, 1., 0., 1.);
	//return clamp(base_with_coverage * cld_coverage, 0, 1);
}

// ----------------------------------------------------------------------------
// Volumetrics
// ----------------------------------------------------------------------------
#define cld_march_steps		(50)
#define cld_thick			(100.)
#define illum_march_steps	(6)

#define sigma_absobtion		(0.)
#define sigma_scattering	(.15)

struct volumetric_sampler_t {
	vec3 origin; // start of ray
	vec3 pos; // current pos of acccumulation ray
	float height;
	float transmittance;
	vec3 radiance;
	float alpha;
};

volumetric_sampler_t define_volume(
	_in(vec3) origin
){
	volumetric_sampler_t v = _begin(volumetric_sampler_t)
		origin,
		origin,
		0.,
		1.,
		vec3(0, 0, 0),
		0.
	_end;
	return v;
}

float illuminate_volumetric(
	_in(vec3) origin,
	_in(float) height,
	_in(vec3) V,
	_in(vec3) L
){
#if 0
	float luminance = exp(height) / 2.;
#else

	const float dt = cld_thick / float(cld_march_steps);
	volumetric_sampler_t vol = define_volume(origin);
	vol.pos += L * dt; // don't sample just where the main raymarcher is

#ifdef HLSL
	[unroll(illum_march_steps)]
#endif
	for (int i = 0; i < illum_march_steps; i++) {
		float density = density_func(vol.pos, 0.);

		vol.transmittance *= exp(- density * sigma_scattering * dt);

		vol.pos += L * dt;
	}

	float luminance = vol.transmittance;
#endif

#if 0
	return luminance;
#else
	return luminance * henyey_greenstein_phase_func(clamp(dot(L, V), 0., 1.));
#endif
}

void integrate_volumetric(
	_inout(volumetric_sampler_t) vol,
	_in(vec3) V,
	_in(vec3) L,
	_in(float) density,
	_in(float) dt
){
	// change in transmittance (follows Beer-Lambert law)
	float T_i = exp(-density * sigma_scattering * dt);
	// Update accumulated transmittance
	vol.transmittance *= T_i;

	// integrate output radiance (here essentially color)
	vol.radiance += 
		(density * sigma_scattering) * 
		illuminate_volumetric(vol.pos, vol.height, V, L) *
		vol.transmittance * 
		8. *
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

#ifdef SPHERE
	hit_t hit = no_hit;
	intersect_sphere(eye, atmosphere, hit);

	vec3 projection = eye.direction;
	vec3 origin = hit.origin;
#else
	vec3 projection = eye.direction / eye.direction.y;
	vec3 origin = eye.origin + projection * 150.;
#endif

	vec3 iter = projection * march_step;
	volumetric_sampler_t cloud = define_volume(origin);

#ifdef HLSL
	[fastopt] [loop]
#endif
	for (int i = 0; i < cld_march_steps; i++) {
		cloud.height = float(i) / float(cld_march_steps);
		
		float density = density_func(cloud.pos, cloud.height);

		integrate_volumetric(
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