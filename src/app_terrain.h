#include "def.h"

struct sphere_t {
	vec3 origin;
	float radius;
	int material;
};

bool intersect_sphere(_in(ray_t) ray, _in(sphere_t) sphere, _inout(float) t0, _inout(float) t1)
{
	vec3 rc = sphere.origin - ray.origin;
	float radius2 = sphere.radius * sphere.radius;
	float tca = dot(rc, ray.direction);
	//	if (tca < 0.) return false;
	float d2 = dot(rc, rc) - tca * tca;
	if (d2 > radius2) return false;
	float thc = sqrt(radius2 - d2);
	t0 = tca - thc;
	t1 = tca + thc;

	if (t0 > t1) {
		float temp = t0;
		t0 = t1;
		t1 = temp;
	}

	if (t0 < 0.) {
		t0 = t1;
		if (t0 < 0.) return false;
	}

	return true;
}

//
// Main
//

const vec3 betaR = vec3(5.5e-6, 13.0e-6, 22.4e-6); // Rayleigh scattering coefficients at sea level (m)
const vec3 betaM = vec3(21e-6); // Mie scattering coefficients at sea level (m)
const float Rh = 7994.0; // Rayleigh scale height (m)
const float Mh = 1200.0; // Mie scale height (m)
const float earth_radius = 6360e3; // (m)
const float atmosphere_radius = 6420e3; // (m)
const vec3 sun_dir = vec3(0, 1, 0);
const float sun_power = 20.0;
const float mean_cos = 0.76; // defines if the light is mainly scattered along the forward or backwards direction

const int air = 1;
const sphere_t atmosphere = sphere_t _begin
	vec3(0), atmosphere_radius, air
_end;

vec3 get_incident_light(_in(ray_t) ray)
{
	float t0, t1;
	if (!intersect_sphere(ray, atmosphere, t0, t1)) return vec3(0);
	return vec3(1, 0, 0);
}

void mainImage(_out(vec4) fragColor, _in(vec2) fragCoord)
{
	vec2 aspect_ratio = vec2(iResolution.x / iResolution.y, 1);
	float fov = tan(radians(45.0));
	vec2 point_ndc = fragCoord.xy / iResolution.xy;
	vec3 point_cam = vec3((2.0 * point_ndc - 1.0) * aspect_ratio * fov, -1.0);

	vec3 col = vec3(0);

	// sky dome angles
	vec3 p = point_cam;
	float z2 = p.x * p.x + p.y * p.y;
	if (z2 <= 1.0) {
		float phi = atan(p.y, p.x); // this is actually atan2 from C
		float theta = acos(1.0 - z2);

		vec3 dir = vec3(
			sin(theta) * cos(phi),
			cos(theta),
			sin(theta) * sin(phi));

		ray_t ray = ray_t _begin
			vec3(0, earth_radius + 1., 0),
			dir
			_end;
		col =
			//abs(dir);
			get_incident_light(ray);
	}

	fragColor = vec4(col, 1);
}
