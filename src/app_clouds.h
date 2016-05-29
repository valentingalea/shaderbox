#include "def.h"
#include "util.h"
#include "volumetric.h"

#include "../lib/ashima-noise/src/common.glsl"
#include "../lib/ashima-noise/src/classicnoise3d.glsl"
#include "../lib/ashima-noise/src/noise3d.glsl"
#include "../lib/ashima-noise/src/cellular3d.glsl"
#define noise(x) (snoise(x))

#include "fbm.h"
DECL_FBM_FUNC(fbm_simple_but_nice, 5, abs(snoise(p)))
DECL_FBM_FUNC(fbm_low_freq_perlin, 4, cnoise(p))
DECL_FBM_FUNC(fbm_low_freq_worley, 3, cellular(p).r)

//#define TEX
#ifdef TEX
Texture3D u_tex_noise : register(t1);
SamplerState u_sampler0 : register(s0);
#endif

float fbm_final(_in(vec3) pos)
{
#ifdef TEX
	return u_tex_noise.SampleLevel(u_sampler0, pos, 0).r;
#endif

	//return fbm_simple_but_nice(pos * 2.03, 2.64, .5, .5);

	float p = fbm_low_freq_perlin(pos * 4., 2., .5, .5);
	float w = 1. - fbm_low_freq_worley(pos * 4., 4., .5, .5);
	float n = remap(p, -w, 1., 0., 1.);
	return n;
}

// ----------------------------------------------------------------------------
// Clouds
// ----------------------------------------------------------------------------
#define cld_march_steps  (50)
#define cld_coverage     (.3475)
#define cld_thick        (100.)
#define cld_absorb_coeff (1.)
#define cld_wind_dir     vec3(0, 0, -u_time * .2)
#define cld_sun_dir      normalize(vec3(0, abs(sin(u_time * .3)), -1))

//#define SPHERE
#ifdef SPHERE
_constant(sphere_t) atmosphere = _begin(sphere_t)
	vec3(0, -475, 0), 500., 0
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
	look_at = mul(rotate_around_y(u_mouse.x), vec3(0, 0, -1));
}

void setup_scene()
{
}

vec3 render_sky_color(
	_in(vec3) eye_dir
){
	_constant(vec3) sun_color = vec3(1., .7, .55);
	float sun_amount = max(dot(eye_dir, cld_sun_dir), 0.);

	vec3 sky = mix(vec3(.0, .1, .4), vec3(.3, .6, .8), 1.0 - eye_dir.y);
	sky += sun_color * min(pow(sun_amount, 1500.0) * 5.0, 1.0);
	sky += sun_color * min(pow(sun_amount, 10.0) * .6, 1.0);

	return abs(sky);
}

float density_func(
	_in(vec3) pos,
	_in(float) h
){
	vec3 p = pos * cld_noise_factor + cld_wind_dir;
	float base = fbm_final(p);

// my old method
//	return base * smoothstep (cld_coverage, cld_coverage + .035, base);

// book equiv method
//	return smoothstep(cld_coverage, 1., base);

// GPU Pro 7
	float base_with_coverage = remap(base, cld_coverage, 1., 0., 1.);
	return clamp(base_with_coverage * cld_coverage, 0, 1);
}

float illuminate_volume(
	_inout(volume_sampler_t) cloud,
	_in(vec3) V,
	_in(vec3) L
){
	return 1.; // exp(cloud.height) / 2.02;
}

void intersect_sphere(
	_in(ray_t) ray,
	_in(sphere_t) sphere,
	_inout(hit_t) hit
) {
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

vec4 render_clouds(
	_in(ray_t) eye
){
	const int steps = cld_march_steps;
	const float march_step = cld_thick / float(steps);

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
	volume_sampler_t cloud = begin_volume(origin, cld_absorb_coeff);

#ifdef HLSL
	[fastopt] [loop]
#endif
	for (int i = 0; i < steps; i++) {
		cloud.height = float(i) / float(steps);
		float dens = density_func(cloud.pos, cloud.height);

		integrate_volume(
			cloud,
			eye.direction, cld_sun_dir,
			dens, march_step);

		cloud.pos += iter;

		if (cloud.alpha > .999) break;
	}

	float cutoff = dot(eye.direction, vec3(0, 1, 0));
	return vec4(cloud.C, cloud.alpha * smoothstep(.0, .2, cutoff));
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

#define FOV 1.//tan(radians(30.))
#include "main.h"