// ----------------------------------------------------------------------------
// Rayleigh and Mie scattering atmosphere system
//
// implementation of the techniques described here:
// http://www.scratchapixel.com/old/lessons/3d-advanced-lessons/simulating-the-colors-of-the-sky/atmospheric-scattering/
// ----------------------------------------------------------------------------

#include "def.h"
#include "util.h"
#include "intersect.h"

#define hg_g (.76)
#include "volumetric.h"

bool isect_sphere(_in(ray_t) ray, _in(sphere_t) sphere, _inout(float) t0, _inout(float) t1)
{
	vec3 rc = sphere.origin - ray.origin;
	float radius2 = sphere.radius * sphere.radius;
	float tca = dot(rc, ray.direction);
	float d2 = dot(rc, rc) - tca * tca;
	float thc = sqrt(radius2 - d2);
	t0 = tca - thc;
	t1 = tca + thc;

	return d2 < radius2;
}

// scattering coefficients at sea level (m)
_constant(vec3) betaR = vec3(5.5e-6, 13.0e-6, 22.4e-6); // Rayleigh 
_constant(vec3) betaM = vec3(21e-6, 21e-6, 21e-6); // Mie

// scale height (m)
// thickness of the atmosphere if its density were uniform
_constant(float) hR = 7994.0; // Rayleigh
_constant(float) hM = 1200.0; // Mie

_constant(float) earth_radius = 6360e3; // (m)
_constant(float) atmosphere_radius = 6420e3; // (m)

_mutable(vec3) sun_dir = vec3(0, 1, 0);
_constant(float) sun_power = 20.0;

_constant(sphere_t) atmosphere = _begin(sphere_t)
	vec3(0, 0, 0), atmosphere_radius, 0
_end;

_constant(int) num_samples = 16;
_constant(int) num_samples_light = 8;

bool get_sun_light(
	_in(ray_t) ray,
	_inout(float) optical_depthR,
	_inout(float) optical_depthM
){
	float t0, t1;
	isect_sphere(ray, atmosphere, t0, t1);

	float march_pos = 0.;
	float march_step = t1 / float(num_samples_light);

	for (int i = 0; i < num_samples_light; i++) {
		vec3 sample =
			ray.origin +
			ray.direction * (march_pos + 0.5 * march_step);
		float height = length(sample) - earth_radius;
		if (height < 0.)
			return false;

		optical_depthR += exp(-height / hR) * march_step;
		optical_depthM += exp(-height / hM) * march_step;

		march_pos += march_step;
	}

	return true;
}

vec3 get_incident_light(_in(ray_t) ray)
{
	// "pierce" the atmosphere with the viewing ray
	float t0, t1;
#ifdef HLSL
	[flatten]
#endif
	if (!isect_sphere(
		ray, atmosphere, t0, t1)) {
		return vec3(0., 0., 0.);
	}

	float march_step = t1 / float(num_samples);

	// cosine of angle between view and light directions
	float mu = dot(ray.direction, sun_dir);

	// Rayleigh and Mie phase functions
	// A black box indicating how light is interacting with the material
	// Similar to BRDF except
	// * it usually considers a single angle
	//   (the phase angle between 2 directions)
	// * integrates to 1 over the entire sphere of directions
	float phaseR = rayleigh_phase_func(mu);
	float phaseM =
#if 1
		henyey_greenstein_phase_func(mu);
#else
		schlick_phase_func(mu);
#endif

	// optical depth (or "average density")
	// represents the accumulated extinction coefficients
	// along the path, multiplied by the length of that path
	float optical_depthR = 0.;
	float optical_depthM = 0.;

	vec3 sumR = vec3(0, 0, 0);
	vec3 sumM = vec3(0, 0, 0);
	float march_pos = 0.;

	for (int i = 0; i < num_samples; i++) {
		vec3 sample =
			ray.origin +
			ray.direction * (march_pos + 0.5 * march_step);
		float height = length(sample) - earth_radius;

		// integrate the height scale
		float hr = exp(-height / hR) * march_step;
		float hm = exp(-height / hM) * march_step;
		optical_depthR += hr;
		optical_depthM += hm;

		// gather the sunlight
		ray_t light_ray = _begin(ray_t)
			sample,
			sun_dir
		_end;
		float optical_depth_lightR = 0.;
		float optical_depth_lightM = 0.;
		bool overground = get_sun_light(
			light_ray,
			optical_depth_lightR,
			optical_depth_lightM);

		if (overground) {
			vec3 tau =
				betaR * (optical_depthR + optical_depth_lightR) +
				betaM * 1.1 * (optical_depthM + optical_depth_lightM);
			vec3 attenuation = exp(-tau);

			sumR += hr * attenuation;
			sumM += hm * attenuation;
		}

		march_pos += march_step;
	}

	return
		sun_power *
		(sumR * phaseR * betaR +
		sumM * phaseM * betaM);
}

#define FROM_SPACE 1

void setup_camera(
	_inout(vec3) eye,
	_inout(vec3) look_at
){
#ifdef FROM_SPACE
	eye = vec3(0, 0, 0);
	look_at = vec3(0, 1, 0);
#else
	eye = vec3(0, earth_radius + 1., 0);
	look_at = vec3(0, earth_radius + 1.5, -1);
#endif
}

void setup_scene()
{
	mat3 rot = rotate_around_x(-abs(sin(u_time / 2.)) * 90.);
	sun_dir = mul(sun_dir, rot);
}

vec3 render(
	_in(ray_t) eye,
	_in(vec3) point_cam
){
	vec3 col = vec3(0, 0, 0);


#ifdef FROM_SPACE
#ifdef HLSL
#define atan(y, x) atan2(x, y)
#endif
	// sky dome angles
	vec3 p = point_cam;
	float z2 = p.x * p.x + p.y * p.y;
	float phi = atan(p.y, p.x);
	float theta = acos(1.0 - z2);
	vec3 dir = vec3(
		sin(theta) * cos(phi),
		cos(theta),
		sin(theta) * sin(phi));

	ray_t ray = _begin(ray_t)
		vec3(0, earth_radius + 1., 0),
		dir
	_end;
	
	col = get_incident_light(ray);
#else
	plane_t terrain = _begin(plane_t)
		vec3 (0, -1, 0),
		earth_radius,
		0
	_end;
	
	hit_t hit = no_hit;
	intersect_plane (eye, terrain, hit);
	
	if (hit.t > max_dist) {
		col = get_incident_light(eye);
	} else {
		col = vec3 (.33, .33, .33);
	}
#endif

	return col;
}

#define FOV 1. // 45 degrees
#include "main.h"
