#include "def.h"
#include "util.h"
#include "intersect.h"
#include "volumetric.h"

#define cld_march_steps (25)
#define cld_coverage (.5)
#define cld_thick (25.)
#define cld_absorb_coeff (1.03)
#define cld_wind_dir vec3(0, 0, -u_time * .2)
#define cld_sun_dir normalize(vec3(0, abs(sin(u_time * .3)), -1))

#if 1
#include "noise_iq.h"
#define noise(x) noise_iq(x)
#else
#include "noise_worley.h"
#define noise(x) (1. - noise_w(x).r)
#endif

#include "fbm.h"
DECL_FBM_FUNC(fbm_clouds, 4, noise(p))

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
	_in(vec3) pos
){
	vec3 p = pos * .0212242 + cld_wind_dir;
	float dens = fbm_clouds(p, 2.76434, .5, .5);
	
	dens *= smoothstep (cld_coverage, cld_coverage + .035, dens);

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
	const float cutoff = dot(eye.direction, vec3(0, 1, 0));

	volume_sampler_t cloud = begin_volume(eye.origin, cld_absorb_coeff);
	cloud.pos = cloud.origin + projection * 100.;
	vec3 iter = projection / march_step;

	for (int i = 0; i < steps; i++) {
		float dens = density_func(cloud.pos);

		cloud.height = float(i) / float(steps);
		integrate_volume(
			cloud,
			eye.direction, cld_sun_dir,
			dens, march_step);

		cloud.pos += iter;
	}

	return vec4(cloud.C, cloud.alpha * smoothstep(.0, .2, cutoff));
}

void setup_camera(
	_inout(vec3) eye,
	_inout(vec3) look_at
){
	eye = vec3(0, 1., 0);
	look_at = vec3(0, 1.6, -1);
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