#include "main_pre.h"
#define SCREEN_WIDTH 300
#define SCREEN_HEIGHT 300
/// GLSL begin //////////////////////////////////////////////////////////////////
#ifdef __cplusplus
#define _in(T) const T &
#define _inout(T) T &
#define _out(T) T &
#define _begin {
#define _end }
#define _rvalue_ref(T) T
#define _move(T) T
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
	float metallic;
	float roughness;
	float ior; // index of refraction
	float reflectivity;
	float translucency;
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

#define num_planes 6
plane_t planes[num_planes];

#define num_spheres 3
sphere_t spheres[num_spheres];

#define num_lights 1
point_light_t lights[num_lights];
vec3 ambient_light = vec3 (.01, .01, .01);

#define num_materials 10
#define mat_invalid -1
#define mat_debug 0
material_t materials[num_materials];

vec3 eye;

// the API
void intersect_sphere (_in(ray_t) ray, _in(sphere_t) sphere, _inout(hit_t) hit);
void intersect_plane (_in(ray_t) ray, _in(plane_t) p, _inout(hit_t) hit);
float fresnel_factor (_in(float) n1, _in(float) n2, _in(float) VdotH);
mat3 rotate_around_y (_in(float) angle_degrees);
mat3 rotate_around_x (_in(float) angle_degrees);
vec3 corect_gamma (_in(vec3) color);
_rvalue_ref(material_t) get_material (_in(int) index);

// GLSL/HLSL utilities ported to C++
float saturate(_in(float) value) { return clamp(value, 0., 1.); }
#ifdef __cplusplus
vec3 faceforward(_in(vec3) N, _in(vec3) I, _in(vec3) Nref);
vec3 reflect(_in(vec3) incident, _in(vec3) normal);
vec3 refract(_in(vec3) incident, _in(vec3) normal, _in(float) n);
#endif

void setup_scene ()
{
	materials [mat_debug] = material_t _begin vec3 (1., 1., 1.), 0., 0., 1., 0., 0. _end;	

//
// Cornell box
//
#define cb_mat_white 1
#define cb_mat_red 2
#define cb_mat_blue 3
#define cb_mat_reflect 4
#define cb_mat_refract 5
	materials[cb_mat_white] = material_t _begin vec3(0.7913), .0, .5, 1., 0., 0. _end;
	materials[cb_mat_red] = material_t _begin vec3(0.6795, 0.0612, 0.0529), 0., .5, 1., 0., 0. _end;
	materials[cb_mat_blue] = material_t _begin vec3(0.1878, 0.1274, 0.4287), 0., .5, 1., 0., 0. _end;
	materials[cb_mat_reflect] = material_t _begin vec3(.5), .0, .0, 1.0, 1., 0. _end;
	materials[cb_mat_refract] = material_t _begin vec3(.5), .0, .0, 1.5, 0., 1. _end;

#define cb_plane_ground 0
#define cb_plane_behind 1
#define cb_plane_front 2
#define cb_plane_ceiling 3
#define cb_plane_left 4
#define cb_plane_right 5
#define cb_plane_dist 2.
	planes[cb_plane_ground] = plane_t _begin vec3(0, -1, 0), 0., cb_mat_white _end;
	planes[cb_plane_ceiling] = plane_t _begin vec3(0, 1, 0), 2. * cb_plane_dist, cb_mat_white _end;
	planes[cb_plane_behind] = plane_t _begin vec3(0, 0, -1), -cb_plane_dist, cb_mat_white _end;
	planes[cb_plane_front] = plane_t _begin vec3(0, 0, 1), cb_plane_dist, cb_mat_white _end;
	planes[cb_plane_left] = plane_t _begin vec3(1, 0, 0), cb_plane_dist, cb_mat_red _end;
	planes[cb_plane_right] = plane_t _begin vec3(-1, 0, 0), -cb_plane_dist, cb_mat_blue _end;

#define cb_sphere_light 0
#define cb_sphere_left 1
#define cb_sphere_right 2
	spheres[cb_sphere_light] = sphere_t _begin vec3(0, 2.5 * cb_plane_dist + 0.4, 0), 1.5, mat_debug _end;
	spheres[cb_sphere_left] = sphere_t _begin vec3(0.75, 1, -0.75), 1., cb_mat_reflect _end;
	spheres[cb_sphere_right] = sphere_t _begin vec3(-0.75, 0.75, 0.75), 0.75, cb_mat_refract _end;

	lights[0] = point_light_t _begin vec3(0, 2. * cb_plane_dist - 0.2, 0), vec3 (1., 1., 1.) _end;
}

vec3 background(_in(ray_t) ray)
{
#if 0
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
	float fresnel_term = fresnel_factor (1., mat.ior, VdotH);

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

hit_t raytrace_iteration (_in(ray_t) ray, _in(int) mat_ignored)
{
	hit_t hit = no_hit;

	for (int i = 0; i < num_planes; ++i) {
		intersect_plane(ray, planes [i], hit);
	}

	for (int i = 0; i < num_spheres; ++i) {
		if (spheres [i].material != mat_invalid
		&& spheres [i].material != mat_ignored) {
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
	hit_t hit = raytrace_iteration (ray, mat_invalid);

	if (hit.t >= max_dist) {
		return background (ray);
	}

#ifdef __cplusplus
#if 1
	_rvalue_ref(material_t) mat = get_material(hit.material_id);

	if ((mat.reflectivity > 0. || mat.translucency > 0.) && depth < MAX_DEPTH) {
		vec3 color = vec3(0);
		
		//TODO: Proper Fresnel
		//float facing_ratio = -dot (ray.direction, hit.normal);
        //float fresnel_effect = mix(pow(1 - facing_ratio, 2), 1, 0.1);
        //return vec3 (fresnel_effect);
        float kr = 1.; //fresnel_effect;
        float kt = 1.; //- kr;

	// reflection
		if (kr * mat.reflectivity > 0.) {
			vec3 refl_dir = reflect(ray.direction, hit.normal);
			ray_t refl_ray = ray_t _begin
				hit.origin + refl_dir * BIAS,
				refl_dir
			_end;
			color += kr * mat.reflectivity *
			raytrace_all(refl_ray, depth + 1);
		}

	// refraction (or transmission)
		if (kt * mat.translucency > 0.) {
			bool outside = dot (ray.direction, hit.normal) < 0;
			float eta = 1. / mat.ior;
			vec3 trans_dir;
			if (outside) {
				trans_dir = normalize (refract (ray.direction, hit.normal, eta));
			}
			else {
				eta = 1. / eta;
				trans_dir = normalize (refract (ray.direction, -hit.normal, eta));
			}
			ray_t trans_ray = ray_t _begin
				hit.origin + trans_dir * BIAS,
				trans_dir
			_end;
			color += kt * mat.translucency *
			raytrace_all (trans_ray, depth + 1);
		}

		return color;
	}
#endif
#endif

	vec3 color = illuminate (hit);

#if 1 // shadow ray
	vec3 sh_line = lights [0].origin - hit.origin;
	vec3 sh_dir = normalize (sh_line);
	ray_t sh_trace = ray_t _begin
		hit.origin + sh_dir * BIAS,
		sh_dir
	_end;
	hit_t sh_hit = raytrace_iteration (sh_trace, mat_debug);
	if (sh_hit.t < length (sh_line)) {
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
	vec2 mouse = iMouse.x < BIAS ? vec2(0) : iResolution.xy / iMouse.xy;
	mat3 rot_y = rotate_around_y(mouse.x * 10.);
	eye = rot_y * vec3 (0, cb_plane_dist, 2.333 * cb_plane_dist);
	vec3 look_at = vec3(0, cb_plane_dist, 0);

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
			hit.material_param = 1.; // cb; // Disabled for now
			hit.origin = impact;
			hit.normal = faceforward(p.direction, ray.direction, p.direction);
		}
	}
}

float fresnel_factor(_in(float) n1, _in(float) n2, _in(float) VdotH)
{
// using Schlick’s approximation    
	float Rn = (n1 - n2) / (n1 + n2);
	float R0 = Rn * Rn; // reflection coefficient for light incoming parallel to the normal
	float F = 1. - VdotH;
	return R0 - (1. - R0) * (F * F * F * F * F);
}

#ifdef __cplusplus
vec3 faceforward(_in(vec3) N, _in(vec3) I, _in(vec3) Nref)
{
	return dot(Nref, I) < 0 ? N : -N;
}

vec3 reflect(_in(vec3) incident, _in(vec3) normal)
{
	return incident - 2. * dot(normal, incident) * normal;
}

vec3 refract (_in(vec3) incident, _in(vec3) normal, _in(float) n)
{
	float cosi = -dot ( normal, incident);
	float sint2 = n * n * (1. - cosi * cosi);
	if (sint2 > 1.) {
		return reflect (incident, normal); // Total Internal Reflection - TODO: is this ok?
	}
	return n * incident + (n * cosi - sqrt (1. - sint2)) * normal;
}
#endif

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

//
// Notes & documentation
//

//     R       V    N    H      L         L dir to light       
//      ^      ^    ^    ^     ^          V dir to eye
//        .     \   |   /    .            N normal
//          .    \  |  /   .              H half between L and V
//            .   \ | /  .                R reflected
//  n1          .  \|/ .                  O hit point             
// -----------------O----------------     T refracted
//  n2             .                      n1 index of refraction of outgoing medium
//                .                       n2 index of refraction of incoming medium
//               .
//              .
//             .
//           \/_ T
//

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
// http://graphics.stanford.edu/courses/cs148-10-summer/docs/2006--degreve--reflection_refraction.pdf
//    
// indices of refraction    
// gold ---- air ---- ice --- water --- beer --- alumim --- glass --- salt --- PET --- asphalt --- lead --- diamond ---- iron ---
// 0.470    1.000    1.309    1.333     1.345     1.390     1.500    1.516    1.575     1.645      2.010     2.420      2.950
//

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