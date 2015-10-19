// https://upload.wikimedia.org/wikipedia/commons/thumb/5/57/Cloud_types_en.svg/960px-Cloud_types_en.svg.png
// http://oceanservice.noaa.gov/education/yos/resource/JetStream/synoptic/clouds_max.htm#max
// http://www.metoffice.gov.uk/learning/clouds/cloud-spotting-guide
// http://www.srh.noaa.gov/srh/jetstream/clouds/cloudwise/types.html

#include "def.h"
#include "noise_iq.h"
#include "fbm.h"
#include "util.h"

vec3 render_sky_color(_in(ray_t) eye, _in(vec3) sun_dir)
{
	float sun_dot = clamp(dot(eye.direction, sun_dir), 0., 1.);

	// colour scheme taken from https://www.shadertoy.com/view/MlSSR1
	vec3 blue = vec3(0.3, .55, 0.8);
	vec3 red = vec3(0.8, 0.8, 0.6);
	vec3 sky = mix(blue, red, 2.45*pow(sun_dot, 228.));

	return sky;// * (1. - 0.18*eye.direction);
}

vec3 render_clouds(_in(ray_t) eye)
{
	if (eye.direction.y < .1) return vec3(0);

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
	vec3 C = vec3(0); // color
	const float absorption = 1.25;
	float alpha = 0.;

	for (int i = 0; i < steps; i++) {
		// sample point
		vec3 p = pos * .0002242 + vec3(0, 0, u_time * .5);

		// density function
		float dens = fbm(p);
		dens = smoothstep(.526, 1., dens);

		// Beer's law
		float T_i = exp(-absorption * dens * march_step);

		T *= T_i;
		if (T < .01)
			break; //return vec3((float(i) * 4.) / 255., 0, 0);

		C += T * /*light(p) */ /*color*/ dens * march_step;
		alpha += (1. - T_i) * (1. - alpha);

		pos += dir_step;
	}

	C *= alpha;

	// add horizon (hide lower artifact/reflection)
	// linear interp the Y with pow func that 
	// ramps up fast at the end, 0 otherwise
	return C;//mix(C, vec3(0), pow(1. - max(eye.direction.y, 0.), 8.));
}

void mainImage(_out(vec4) fragColor, _in(vec2) fragCoord)
{
	vec2 aspect_ratio = vec2(u_res.x / u_res.y, 1);
	float fov = tan(radians(45.0));
	vec2 point_ndc = fragCoord.xy / u_res.xy;
	vec3 point_cam = vec3((2.0 * point_ndc - 1.0) * aspect_ratio * fov, -1.0);

	vec3 col = vec3(0, 0, 0);

	vec3 sun_dir = vec3(0, 0, -1);
//	mat3 rot = rotate_around_x(-abs(sin(u_time / 2.)) * 90.);
//	sun_dir *= rot;

	vec3 eye = vec3(0, 1., 0);
	vec3 look_at = vec3(0, 1.5, -1);
	ray_t eye_ray = get_primary_ray(point_cam, eye, look_at);

//	col += render_sky_color(eye_ray, sun_dir);
	col += render_clouds(eye_ray);

	fragColor = vec4(col, 1);
}