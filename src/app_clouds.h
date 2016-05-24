#include "def.h"
#include "util.h"
#include "intersect.h"
#include "volumetric.h"

#define cld_march_steps  (50)
#define cld_coverage     (.3125)
#define cld_thick        (50.)
#define cld_dist         (50.)
#define cld_absorb_coeff (1.)
#define cld_wind_dir vec3(0, 0, -u_time * .2)
#define cld_sun_dir normalize(vec3(0, abs(sin(u_time * .3)), -1))
_mutable(float) coverage_map;

#include "noise_iq.h"
#define gnoise(x) noise_iq(x)

#include "../lib/ashima-noise/src/noise3d.glsl"
#define noise(x) snoise(x)

#include "fbm.h"
DECL_FBM_FUNC(fbm_clouds, 5, abs(noise(p)))

vec3 render_sky_color(
	_in(vec3) eye_dir
){
	_constant(vec3) sun_color = vec3(1., .7, .55);
	float sun_amount = max(dot(eye_dir, cld_sun_dir), 0.);

	vec3 sky = mix(vec3(.0, .1, .4), vec3(.3, .6, .8), 1.0 - eye_dir.y);
	sky += sun_color * min(pow(sun_amount, 1500.0) * 5.0, 1.0);
	sky += sun_color * min(pow(sun_amount, 10.0) * .6, 1.0);

	return sky;
}

float density_func(
	_in(vec3) pos,
	_in(float) h
){
	vec3 p = pos / (cld_dist * 10.) + cld_wind_dir;
	float dens = fbm_clouds(p * 2.03, 2.64, .5, .5);
	
	dens *= smoothstep (cld_coverage, cld_coverage + .035, dens);

	//dens *= band(.2, .3, .5 + coverage_map * .5, h);

	return dens;
}

float illuminate_volume(
	_inout(volume_sampler_t) cloud,
	_in(vec3) V,
	_in(vec3) L
){
	return exp(cloud.height) / 1.75;
}

vec4 render_clouds(
	_in(ray_t) eye
){
	const int steps = cld_march_steps;
	const float march_step = cld_thick / float(steps);

	const vec3 projection = eye.direction / eye.direction.y;
	const vec3 iter = projection * march_step;

	const float cutoff = dot(eye.direction, vec3(0, 1, 0));

	volume_sampler_t cloud = begin_volume(
		eye.origin + projection * cld_dist,
		cld_absorb_coeff);

	//coverage_map = gnoise(projection);
	//return vec4(coverage_map, coverage_map, coverage_map, 1);

	for (int i = 0; i < steps; i++) {
		cloud.height = (cloud.pos.y - cloud.origin.y)
			/ cld_thick;
		float dens = density_func(cloud.pos, cloud.height);

		integrate_volume(
			cloud,
			eye.direction, cld_sun_dir,
			dens, march_step);

		cloud.pos += iter;

		if (cloud.alpha > .999) break;
	}

	return vec4(cloud.C, cloud.alpha * smoothstep(.0, .2, cutoff));
}

void setup_camera(
	_inout(vec3) eye,
	_inout(vec3) look_at
){
	eye = vec3(0, 0, 0);
	look_at = vec3(0, .75, -1);
}

void setup_scene()
{
}

vec3 render(
	_in(ray_t) eye_ray,
	_in(vec3) point_cam
){
	vec3 sky = render_sky_color(eye_ray.direction);
	if (dot(eye_ray.direction, vec3(0, 1, 0)) < 0.05) return sky;

	vec4 cld = render_clouds(eye_ray);
	vec3 col = mix(sky, cld.rgb, cld.a);

	return col;
}

#define FOV 1. // 45 degrees
#include "main.h"