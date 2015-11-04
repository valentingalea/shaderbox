// https://upload.wikimedia.org/wikipedia/commons/thumb/5/57/Cloud_types_en.svg/960px-Cloud_types_en.svg.png
// http://oceanservice.noaa.gov/education/yos/resource/JetStream/synoptic/clouds_max.htm#max
// http://www.metoffice.gov.uk/learning/clouds/cloud-spotting-guide
// http://www.srh.noaa.gov/srh/jetstream/clouds/cloudwise/types.html

#include "def.h"
#include "noise_iq.h"
#include "fbm.h"
#include "util.h"
#include "intersect.h"

_mutable(vec3) sun_dir = vec3(0, 0, -1);

vec3 render_sky_color(
	_in(ray_t) eye
){
	float sun_dot = clamp(dot(eye.direction, sun_dir), 0., 1.);

	vec3 blue = vec3(0.25, .55, 0.85);
	vec3 red = vec3(0.9, 0.9, 0.75);
	vec3 sky = mix(blue, red, 2.5 * pow(sun_dot, 228.));

	return sky;
}

_constant(float) absorption = 1.25;

float density(
	_in(vec3) pos,
	_in(vec3) offset
){
	vec3 p = pos * .0002242 + offset;
	float dens = fbm(p);
	return smoothstep(.526, 1., dens);	
}

// TODO: not working correctly
float light(
	_in (vec3) origin
){
	const int steps = 2;
	float march_step = .025;
	
	vec3 pos = origin;
	vec3 dir_step = sun_dir * march_step;
	
	float T = 1.; // transmitance
	
	for (int i = 0; i < steps; i++) {
		float dens = density (pos, vec3(0, 0, 0));

		float T_i = exp(-absorption * dens * march_step);

		T *= T_i;
		if (T < .01)
			break;
			
		pos += dir_step;
	}
	
	return T;
}

vec4 render_clouds(
	_in(ray_t) eye
){
	const int steps = 64;
	const float thickness = 100.;
	float march_step = thickness / float(steps);

	// create vanishing point by projection with Y
	const float sky_height = 200.0;
	float dist = sky_height / eye.direction.y;
	vec3 dir = eye.direction * dist;

	vec3 dir_step = dir * march_step;
	vec3 pos = eye.origin;

	float T = 1.; // transmitance
	vec3 C = vec3(0, 0, 0); // color
	float alpha = 0.;

	for (int i = 0; i < steps; i++) {
		float dens = density (pos, vec3(0, 0, u_time * .5));

		float T_i = exp(-absorption * dens * march_step);

		T *= T_i;
		if (T < .01)
			break; //return vec3((float(i) * 4.) / 255., 0, 0);

		C += T * /*light(pos)*/ /*color*/ dens * march_step;
		alpha += (1. - T_i) * (1. - alpha);

		pos += dir_step;
		if (length(pos) > 1e5)
			break;
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
	vec3 point_cam = vec3((2.0 * point_ndc - 1.0) * aspect_ratio * fov, -1.0);

	vec3 col = vec3(0, 0, 0);

//	mat3 rot = rotate_around_x(-abs(sin(u_time / 2.)) * 90.);
//	sun_dir *= rot;

	vec3 eye = vec3(0, 1., 0);
	vec3 look_at = vec3(0, 1.5, -1);
	ray_t eye_ray = get_primary_ray(point_cam, eye, look_at);

//	eye_ray.direction.yz *= rotate_2d(+u_mouse.y * .13);
//	eye_ray.direction.xz *= rotate_2d(-u_mouse.x * .33);

	plane_t ground = _begin(plane_t)
		vec3(0., -1., 0.), 0., 0
		_end;
	hit_t hit = no_hit;
	intersect_plane(eye_ray, ground, hit);

	if (hit.t < max_dist) {
		float cb = checkboard_pattern(hit.origin.xz, .5);
		col = mix(vec3(.6, .6, .6), vec3(.75, .75, .75), cb);
	} else {
		vec3 sky = render_sky_color(eye_ray);
		vec4 cld = render_clouds(eye_ray);
		col = mix(sky, cld.rgb, cld.a);
	}

	fragColor = vec4(col, 1);
}