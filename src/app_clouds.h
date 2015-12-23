// https://upload.wikimedia.org/wikipedia/commons/thumb/5/57/Cloud_types_en.svg/960px-Cloud_types_en.svg.png
// http://oceanservice.noaa.gov/education/yos/resource/JetStream/synoptic/clouds_max.htm#max
// http://www.metoffice.gov.uk/learning/clouds/cloud-spotting-guide
// http://www.srh.noaa.gov/srh/jetstream/clouds/cloudwise/types.html

#include "def.h"
#include "util.h"
#include "intersect.h"

#include "noise_iq.h"
#include "noise_worley.h"
//#include "../lib/ashima-noise/src/classicnoise3D.glsl"

#define noise(x) noise_iq(x)
//#define noise(x) pnoise(x, vec3(128, 128, 128))
//#define noise(x) (1. - noise_w(x).r)
//#define noise(x) abs( noise_iq(x / 8.) - (1. - (noise_w(x * 2.).r)))
#include "fbm.h"

#ifdef HLSL
Texture3D u_tex_noise : register(t1);
SamplerState u_sampler0 : register(s0);
#endif

float get_noise(_in(vec3) x)
{
#if 0
	return u_tex_noise.Sample(u_sampler0, x);
#else
	return fbm(x);
#endif
}

_constant(vec3) sun_color = vec3(1., .7, .55);
_mutable(vec3) sun_dir = normalize(vec3(0, 0.1, -1));

_mutable(vec3) wind_dir = vec3(0, 0, -u_time * .5);
_constant(float) absorption = .30725;

_constant(sphere_t) atmosphere = _begin(sphere_t)
	vec3(0, -450, 0), 500., 0
_end;
_constant(sphere_t) atmosphere_2 = _begin(sphere_t)
	atmosphere.origin, atmosphere.radius + 50., 0
_end;
_constant(plane_t) ground = _begin(plane_t)
	vec3(0., -1., 0.), 0., 1
_end;

vec3 render_sky_color(
	_in(ray_t) eye
){
	vec3 rd = eye.direction;
	float sun_amount = max(dot(rd, sun_dir), 0.0);

	vec3  sky = mix(vec3(.0, .1, .4), vec3(.3, .6, .8), 1.0 - rd.y);
	sky = sky + sun_color * min(pow(sun_amount, 1500.0) * 5.0, 1.0);
	sky = sky + sun_color * min(pow(sun_amount, 10.0) * .6, 1.0);

	return sky;
}

float density(
	_in(vec3) pos,
	_in(vec3) offset,
	_in(float) t
){
	// signal
	vec3 p = pos * .0212242 +offset;
	float dens = get_noise(p);
	
	//dens = band (.1, .3, .6, dens);
	//dens *= step(.5, dens);
	dens *= smoothstep (.55, .6, dens);
	//dens -= .5;

	return clamp(dens, 0., 1.);	
}

float light(
	_in(vec3) origin
){
	const int steps = 8;
	float march_step = 1.;

	vec3 pos = origin;
	vec3 dir_step = sun_dir * march_step;
	float T = 1.; // transmitance

	for (int i = 0; i < steps; i++) {
		float dens = density(pos, wind_dir, 0.);

		float T_i = exp(-absorption * dens * march_step);
		T *= T_i;
		//if (T < .01) break;

		pos += dir_step;
	}

	return T;
}

vec4 render_clouds(
	_in(ray_t) eye
){
	hit_t hit = no_hit;
	intersect_sphere(eye, atmosphere, hit);

	hit_t hit_2 = no_hit;
	intersect_sphere(eye, atmosphere_2, hit_2);

	const float thickness = 25.; // length(hit_2.origin - hit.origin);
	//const float r = 1. - ((atmosphere_2.radius - atmosphere.radius) / thickness);
	const int steps = 50; // +int(32. * r);
	float march_step = thickness / float(steps);

	vec3 dir_step = eye.direction /* eye.direction.y */ * march_step;
	vec3 pos = //eye.origin + eye.direction * 100.; 
		hit.origin;

	float T = 1.; // transmitance
	vec3 C = vec3(0, 0, 0); // color
	float alpha = 0.;

	for (int i = 0; i < steps; i++) {
		float h = float(i) / float(steps);
		float dens = density (pos, wind_dir, h);

		float T_i = exp(-absorption * dens * march_step);
		T *= T_i;
		//if (T < .01) break;

		C += T * light(pos) * /*color*/ dens * march_step;
		alpha += (1. - T_i) * (1. - alpha);

		pos += dir_step;
		if (length(pos) > 1e3) break;
	}

	return vec4(C, alpha);
}

void mainImage(
	_out(vec4) fragColor,
	_in(vec2) fragCoord
){
	vec2 aspect_ratio = vec2(u_res.x / u_res.y, 1);
	float fov = tan(radians(45.0));
	vec2 point_ndc = fragCoord.xy / u_res.xy;
#ifdef HLSL
	point_ndc.y = 1. - point_ndc.y;
#endif
	vec3 point_cam = vec3((2.0 * point_ndc - 1.0) * aspect_ratio * fov, -1.0);

#if 0
	float n = saturate(get_noise(point_cam));
	fragColor = vec4(vec3(n, n, n), 1);
	return;
#endif

	vec3 col = vec3(0, 0, 0);

	//mat3 rot = rotate_around_x(abs(sin(u_time / 2.)) * 45.);
	//sun_dir = mul(rot, sun_dir);

	vec3 eye = vec3(0, 1., 0);
	vec3 look_at = vec3(0, 1.6, -1);
	ray_t eye_ray = get_primary_ray(point_cam, eye, look_at);

	eye_ray.direction.yz = mul(rotate_2d(+u_mouse.y * .13), eye_ray.direction.yz);
	eye_ray.direction.xz = mul(rotate_2d(-u_mouse.x * .33), eye_ray.direction.xz);

	hit_t hit = no_hit;
	intersect_plane(eye_ray, ground, hit);

	if (hit.material_id == 1) {
		float cb = checkboard_pattern(hit.origin.xz, .5);
		col = mix(vec3(.6, .6, .6), vec3(.75, .75, .75), cb);
	} else {
#if 1
		vec3 sky = render_sky_color(eye_ray);
		vec4 cld = render_clouds(eye_ray);
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