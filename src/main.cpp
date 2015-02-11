#include "main_pre.h"
#define SCREEN_WIDTH 100
#define SCREEN_HEIGHT 100
/// GLSL begin //////////////////////////////////////////////////////////////////
#ifdef __cplusplus
#define _in(T) const T &
#define _inout(T) T &
#define _out(T) T &
#define _begin {
#define _end }
#define _rvalue_ref(T) T &&
#define _move(T) std::move(T)
#else
#define _in(T) const in T
#define _inout(T) inout T
#define _out(T) out T
#define _begin (
#define _end )
#define _rvalue_ref(T) T
#define _move(T) T
#endif

struct ray_t {
	vec3 origin;
	vec3 direction;
};

struct material_t {
	vec3 base_color;
	float roughness;
	float metallic;
	float refraction_index;
};

struct sphere_t {
	vec3 origin;
	float radius;
	int material;
};

struct plane_t {
	vec3 direction;
	float distance;
	int material;
};

struct point_light_t {
	vec3 origin;
	vec3 color;
};

struct hit_t {
	float t;
	int material_id;
	float material_param; // used by dynamic materials to store temp data
	vec3 normal;
	vec3 origin;
};

#define max_dist 999.0
hit_t no_hit = hit_t _begin
	(max_dist + 1.0), -1, 1., vec3 (0), vec3 (0)
_end;

#define num_planes 1
#define ground 0
plane_t planes[num_planes];

#define num_spheres 3
sphere_t spheres[num_spheres];

#define num_lights 1
point_light_t lights[num_lights];
vec3 ambient_light = vec3 (.01, .01, .01);

#define num_materials 3
#define mat_invalid -1
#define mat_debug 0
#define mat_ground 1
#define mat_plastic 2
material_t materials[num_materials];

vec3 eye;

// the API
void intersect_sphere (_in(ray_t) ray, _in(sphere_t) sphere, _inout(hit_t) hit);
void intersect_plane (_in(ray_t) ray, _in(plane_t) p, _inout(hit_t) hit);
float get_fresnel_term (_in(float) n1, _in(float) n2, _in(float) VdotH);
vec3 refract (_in(vec3) incident, _in(vec3) normal, _in(float) n1, _in(float) n2);
mat3 rotate_around_y (_in(float) angle_degrees);
mat3 rotate_around_x (_in(float) angle_degrees);
vec3 corect_gamma (_in(vec3) color);
_rvalue_ref(material_t) get_material (_in(int) index);
float saturate(_in(float) value) { return clamp(value, 0., 1.); }

void setup_scene ()
{
	materials [mat_debug] = material_t _begin
		vec3 (1., 1., 1.), // base color
		1.000, // roughness
		1.000, // metallic
		1.000  // refraction index
	_end;	
	materials [mat_ground] = material_t _begin
		vec3 (.9, .6, .1), // base color
		0.100, // roughness
		1.000, // metallic
		1.300  // refraction index
	_end;
	materials [mat_plastic] = material_t _begin
		vec3 (0., 1., 0.), // base color
		0.100, // roughness
		1.000, // metallic
		1.300  // refraction index 
	_end;

	planes[ground] = plane_t _begin
		vec3 (0, -1, 0), // TODO: why does the normal need to be facing down???
		0.0,
		mat_ground
	_end;

	mat3 rot = rotate_around_y(-50. * iGlobalTime);
	vec3 L = rot * vec3 (0, 1, 2);

	spheres [0] = sphere_t _begin
		vec3 (2, 1, 0),
		1.0,
		mat_invalid
	_end;
	spheres [1] = sphere_t _begin
		vec3 (0, 2, 0),
		1.0,
		mat_plastic
	_end;
	spheres [2] = sphere_t _begin
		L,
		0.2,
		mat_debug
	_end;

	lights [0] = point_light_t _begin
		L,
		vec3 (1., 1., 1.)
	_end;
}

vec3 setup_background(_in(ray_t) ray)
{
#if 1
	vec3 sun_dir = normalize(vec3(0, 0.5 * sin (iGlobalTime/2.), -1)); // sin(iGlobalTime), 1, cos(iGlobalTime)));
	vec3 sun_color = vec3 (2, 2, 0);
	float sun_grad = max(0., dot(ray.direction, sun_dir));
	float sky_grad = dot(ray.direction, vec3 (0, 1, 0));
	return
		pow(sun_grad, 250.) * sun_color +
		pow (sky_grad, 2.) * vec3 (.3, .3, 3.);
#else
	return vec3 (.1, .1, .7);
#endif
}

//     R       V    N    H      L         L dir to light       
//      ^      ^    ^    ^     ^          V dir to eye
//        .     \   |   /    .            N normal
//          .    \  |  /   .              H half between L and V
//            .   \ | /  .                R reflected
//              .  \|/ .                  O hit point             
// -----------------O----------------
//               surface
vec3 illum_point_light_blinn_phong(
    _in(vec3) V, _in(point_light_t) light, _in(hit_t) hit, _in(material_t) mat)
{
	vec3 L = normalize (light.origin - hit.origin);

	vec3 diffuse = max(0., dot (L, hit.normal)) * (mat.base_color * hit.material_param) * light.color;

    float spec_factor = 50.;
#if 0 // Blinn specular
	vec3 H = normalize(L + V);
	vec3 specular = pow (max (0., dot (H, hit.normal)), spec_factor) * light.color; // * specular color
#else // Phong specular
	vec3 R = reflect(-L, hit.normal);
	vec3 specular = pow (max (0., dot (R, V)), spec_factor) * light.color; // * specular color
#endif

	return diffuse + specular;
}

vec3 illum_point_light_cook_torrance(
	_in(vec3) V, _in(point_light_t) light, _in(hit_t) hit, _in(material_t) mat)
{
	vec3 L = normalize(light.origin - hit.origin);
	vec3 H = normalize(L + V);
	float NdotL = saturate (dot (hit.normal, L));
	float NdotH = saturate (dot (hit.normal, H));
	float NdotV = saturate (dot (hit.normal, V));
	float VdotH = saturate (dot (V, H));

	// geometric term
	float geo_term = min(
		1.,
		min(
			(2. * NdotH * NdotV) / VdotH,
			(2. * NdotH * NdotL) / VdotH
		)
	);

	// roughness term
	// using Beckmann Distribution
	float rough_sq = mat.roughness * mat.roughness;
	float rough_a = 1. / (rough_sq * NdotH * NdotH * NdotH * NdotH);
	float rough_exp = (NdotH * NdotH - 1.) / (rough_sq * NdotH * NdotH);
	float rough_term = rough_a * exp(rough_exp);

	// Fresnel term
	float fresnel_term = get_fresnel_term (1., mat.refraction_index, VdotH);

	float specular = (geo_term * rough_term * fresnel_term) / (NdotV * NdotL);
	return NdotL * (specular + (mat.base_color * hit.material_param));
}

vec3 illuminate (_in(hit_t) hit)
{
	_rvalue_ref(material_t) mat = get_material(hit.material_id);

	// special case for debug stuff - just solid paint it
	if (hit.material_id == mat_debug) {
		return materials [mat_debug].base_color;
	}

	vec3 accum = ambient_light; // really cheap equivalent for indirect light 
	vec3 V = normalize (eye - hit.origin); // view direction

	for (int i = 0; i < num_lights; ++i) {
#if 0
		accum += illum_point_light_blinn_phong (V, lights [i], hit, mat);
#else
		accum += illum_point_light_cook_torrance (V, lights [i], hit, mat);
#endif
	}

	return accum;
}

hit_t raytrace_iteration (_in(ray_t) ray)
{
	hit_t hit = no_hit;

	intersect_plane (ray, planes [ground], hit);

	for (int i = 0; i < num_spheres; ++i) {
		if (spheres [i].material != mat_invalid) {
			intersect_sphere(ray, spheres [i], hit);
		}
	}

	return hit;
}

#define BIAS 1e-4 // small offset to add to ray when retracing to avoid self-intersection
#define PI 3.14159265359
#define MAX_DEPTH 3

vec3 raytrace_all (_in(ray_t) ray, _in(int) depth)
{
	hit_t hit = raytrace_iteration (ray);
	
	if (hit.t >= max_dist) {
		return setup_background (ray);
	}
	
#ifdef __cplusplus
#if 0
// reflection
	if (hit.material_id == mat_plastic && depth < MAX_DEPTH) {
		ray_t refl = ray_t _begin
			hit.origin + BIAS,
			reflect (ray.direction, hit.normal)
		_end;
		
		return raytrace_all (refl, depth + 1);
	}
#endif

#if 1
// refraction
	if (hit.material_id == mat_plastic && depth < MAX_DEPTH) {
		vec3 dir;
		float ior_outside = 1.;
		float ior_inside = 1.333;

		if (dot (ray.direction, hit.normal) < 0) {
			dir = normalize (refract (ray.direction, hit.normal, ior_outside, ior_inside));
			//return abs (dir);
		} else {
			dir = normalize (refract (ray.direction, -hit.normal, ior_inside, ior_outside));
			//return abs (dir);
		}

		ray_t trans = ray_t _begin
			hit.origin + dir * (BIAS),
			dir
		_end;
		
		return raytrace_all (trans, depth + 1);
	}
#endif
#endif

	vec3 color = illuminate (hit);

#if 0 // shadow ray
	ray_t trace = ray_t _begin
		hit.origin + BIAS,
		normalize (lights [0].origin - hit.origin)
	_end;
	hit_t sh = raytrace_iteration (trace);
	if (sh.t < length (lights [0].origin - hit.origin)) {
		color *= 0.1;
	}
#endif

	return color;
}

ray_t get_primary_ray (_in(vec3) cam_local_point, _in(vec3) cam_origin, _in(vec3) cam_look_at)
{
	vec3 fwd = normalize (cam_look_at - cam_origin);
	vec3 up = vec3 (0, 1, 0);
	vec3 right = cross (up, fwd);
	up = cross (fwd, right);

	return ray_t _begin
		cam_origin,
		normalize (
			fwd +
			up * cam_local_point.y +
			right * cam_local_point.x
		)
	_end;
}

void main()
{
	// gl_FragCoord is in raster space [0..resolution]
	// convert to NDC [0..1]
	// add 0.5 to select center of "pixel"
	vec2 point_ndc = (gl_FragCoord.xy + 0.5) / iResolution.xy;

	// assuming screen width is larger than height 
	vec2 aspect_ratio = vec2(iResolution.x / iResolution.y, 1);

	// FOV //TODO: more info about this	
	float fov = tan(radians (30.0));

	// convert to camera space [-1, +1]
	// and z is negative because of cross product and the need for x to be to the right 
	// TODO: draw fancy diagram 
	vec3 point_cam = vec3((2.0 * point_ndc - 1.0) * aspect_ratio * fov, -1.0);

	// TODO: make mouse y rotation work
//	vec2 mouse = 2. * (iMouse.x > 0. ? iResolution.xy / iMouse.xy : vec2(0)) - 1.;
//	mat3 rot_y = rotate_around_y(mouse.x * PI * 10.);
	eye = rotate_around_x (iGlobalTime * 10.) * vec3 (0, 0, 5);
	vec3 look_at = vec3(0, 1, 0);

	ray_t ray = get_primary_ray (point_cam, eye, look_at);

	setup_scene();

	vec3 color = raytrace_all (ray, 0);

	gl_FragColor = vec4 (corect_gamma (color), 1);
}

//
// Intersections
//
void intersect_sphere(_in(ray_t) ray, _in(sphere_t) sphere, _inout(hit_t) hit)
{
#if 1
// geometrical solution
// info: http://www.scratchapixel.com/old/lessons/3d-basic-lessons/lesson-7-intersecting-simple-shapes/ray-sphere-intersection/
	vec3 rc = sphere.origin - ray.origin;
	float radius2 = sphere.radius * sphere.radius;
	float tca = dot(rc, ray.direction);
	if (tca < 0.) return;
	float d2 = dot(rc, rc) - tca * tca;
	if (d2 > radius2) return;
	float thc = sqrt(radius2 - d2);
	float t0 = tca - thc;
	float t1 = tca + thc;

	if (t0 < 0.) t0 = t1;
	if (t0 < hit.t) {
#else // TODO: wrong for some reason... t gets weird values at intersection
// analytical solution
// based on combining the
// sphere eq: (P - C)^2 = R^2
// ray eq: P = O + t*D
// into a quadratic eq: ax^2 + bx + c = 0
// which can be solved by "completing the square" http://www.mathsisfun.com/algebra/completing-square.html
// NOTE: be careful about "catastrophic cancellation" http://en.wikipedia.org/wiki/Loss_of_significance
	vec3 rc = ray.origin - sphere.origin;
//	float a = D dot D -- which is 1 because D is normalised
	float b = 2.0 * dot(rc, ray.direction); // 2 * (O - C) dot D
	float c = dot(rc, rc) - sphere.radius * sphere.radius; // (O - C)^2 - R^2
	float discr = b * b - 4.0 * c;
	float t = -b - sqrt(abs(discr)); // use abs to avoid the check for < 0

	if (discr > 0.0 && t > 0.0 && t < hit.t) {
#endif
		vec3 impact = ray.origin + ray.direction * t0;

		hit.t = t0;
		hit.material_id = sphere.material;
		hit.material_param = 1.;
		hit.origin = impact;
		hit.normal = (impact - sphere.origin) / sphere.radius;
	}
}

// Plane is define by normal N and distance to origin P0 (which is on the plane itself)
// a plane eq is: (P - P0) dot N = 0
// which means that any line on the plane is perpendicular to the plane normal
// a ray eq: P = O + t*D
// substitution and solving for t gives:
// t = ((P0 - O) dot N) / (N dot D)
void intersect_plane(_in(ray_t) ray, _in(plane_t) p, _inout(hit_t) hit)
{
	float denom = dot(p.direction, ray.direction);
	if (denom > 1e-6)
	{
		float t = dot(vec3 (p.distance) - ray.origin, p.direction) / denom;
		if (t >= 0.0 && t < hit.t)
		{
			vec3 impact = ray.origin + ray.direction * t;
			vec2 pattern = floor (impact.xz * 0.5);
			float cb = mod (pattern.x + pattern.y, 2.0);

			hit.t = t;
			hit.material_id = p.material;
			hit.material_param = cb;
			hit.origin = impact;
			hit.normal = -p.direction; // TODO: why?
		}
	}
}

// Fresnel effect says that reflection
// increases on a surface with the angle of view
// and the transmission acts the reverse way
// ex: soap bubble - the edges are highly reflected and colorful
// while the middle is transparent
// ex: water of a lake: close to shore you see thru,
// farther away you see only the reflections on the water
//
// how much reflection and how much transmission
// is given by Fresnel equations which can be
// approximated with the Schilck eq
//
// the direction of the refracted (aka transmitted)
// ray is given by the Snell eq given the
// indices of refraction of the mediums
//
// http://www.scratchapixel.com/old/lessons/3d-basic-lessons/lesson-14-interaction-light-matter/optics-reflection-and-refraction/
// http://graphics.stanford.edu/courses/cs148-10-summer/docs/2006--degreve--reflection_refraction.pdf77
//    
// indices of refraction    
// gold ---- air ---- ice --- water --- beer --- alumim --- glass --- salt --- PET --- asphalt --- lead --- diamond ---- iron ---
// 0.470    1.000    1.309    1.333     1.345     1.390     1.500    1.516    1.575     1.645      2.010     2.420      2.950
//
float get_fresnel_term(_in(float) n1, _in(float) n2, _in(float) VdotH)
{
// using Schlick’s approximation    
	float Rn = (n1 - n2) / (n1 + n2);
	float R0 = Rn * Rn; // reflection coefficient for light incoming parallel to the normal
	float F = 1. - VdotH;
	return R0 - (1. - R0) * (F * F * F * F * F);
}

vec3 refract (_in(vec3) incident, _in(vec3) normal, _in(float) n1, _in(float) n2)
{
	float n = n1 / n2;
	float cosi = -dot ( normal, incident);
	float sint2 = n * n * (1. - cosi * cosi);
	if (sint2 > 1.) {
		return reflect (incident, normal); // Total Internal Reflection
	}
	return n * incident + (n * cosi - sqrt (1. - sint2)) * normal;
}

//
// Utils
//
mat3 rotate_around_y(_in(float) angle_degrees)
{
	float angle = radians(angle_degrees);
	float _sin = sin(angle);
	float _cos = cos(angle);
	return mat3(_cos, 0, _sin, 0, 1, 0, -_sin, 0, _cos);
}

mat3 rotate_around_x(_in(float) angle_degrees)
{
	float angle = radians(angle_degrees);
	float _sin = sin(angle);
	float _cos = cos(angle);
	return mat3 (1, 0, 0, 0, _cos, -_sin, 0, _sin, _cos);
}

vec3 corect_gamma(_in(vec3) color)
{
	float gamma = 1.0 / 2.25;
	return vec3 (pow(color.r, gamma), pow(color.g, gamma),pow(color.b, gamma));
}

_rvalue_ref(material_t) get_material(_in(int) index)
{
	_rvalue_ref(material_t) mat = _move(materials[0]);

	for (int i = 0; i < num_materials; ++i) {
		if (i == index) {
			mat = materials[i];
			break;
		}
	}

	return _move(mat);
}
/// GLSL end //////////////////////////////////////////////////////////////////

	// be a dear a clean up
#pragma warning(pop)
#undef main
#undef in
#undef out
#undef inout
#undef uniform
}

// these headers, especially SDL.h & time.h set up names that are in conflict
// with sandbox'es;
// sandbox should be moved to a separate h/cpp pair, but out of laziness...
// including them
// *after* sandbox solves it too

#include <iostream>
#include <sstream>
#include <SDL.h>
#include <SDL_image.h>
#include <time.h>
#include <memory>
#include <functional>
#if OMP_ENABLED
#include <omp.h>
#endif

#include "main_post.h"