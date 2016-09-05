#include "def.h"
#include "util.h"
#include "sdf.h"
#include "material.h"
#include "util_optics.h"
#include "light.h"

// ----------------------------------------------------------------------------
// Vinyl disk animation
// ----------------------------------------------------------------------------

vec3 background(_in(ray_t) ray)
{
	return vec3(.7, .7, .7);
}

#define mat_groove 1
#define mat_dead_wax 2
#define mat_label 3
#define mat_logo 4
#define mat_count (5)

void setup_mat(
	_inout(material_t) mat,
	_in(vec3) diffuse,
	_in(float) metallic,
	_in(float) roughness
){
	mat.base_color = diffuse;
	mat.metallic = metallic;
	mat.roughness = roughness;
	mat.ior = 1.;
	mat.reflectivity = 0.;
	mat.translucency = 0.;
}

void setup_scene()
{
	setup_mat(materials[mat_debug], vec3(1, 1, 1), 0, 0);
	setup_mat(materials[mat_groove], vec3(.1, .1, .1), 0, .13);
	setup_mat(materials[mat_dead_wax], vec3(.2, .2, .2), 0, .2);
	setup_mat(materials[mat_label], vec3(.5, .1, .1), 0, .5);
	setup_mat(materials[mat_logo], vec3(0, 0, .7), 0, .5);
}

void setup_camera(_inout(vec3) eye, _inout(vec3) look_at)
{
	eye = vec3(0, 0, 7);
	look_at = vec3(0, 0, 0);
}

vec2 sdf(_in(vec3) pos)
{
	// NOTE: everything is centered around origin
	// change coord frame by offseting
	// with inverse then doing
	// the opposite before next
	// effectively doing push/pop

	// NOTE: all measurements are in halfs
	// due to the above
	
	vec3 p = mul(pos, rotate_around_x(90));
	p = mul(p, rotate_around_y(u_time * 33.));
	const float thick = .1;
	const float D = .51;
	const float R = 6.;
	
	vec2 groove = vec2(
		sd_y_cylinder(p, 6., thick),
		mat_groove);
	vec2 dead_wax = vec2(
		sd_y_cylinder(p, 2.5, thick + D),
		mat_dead_wax);
	vec2 label = vec2(
		sd_y_cylinder(p, 2., thick + D * 2.),
		mat_label);
	vec2 logo = vec2(
		sd_box(p, vec3(1, thick + D, 1)),
		mat_logo);
	float center_hole =
		sd_y_cylinder(p, .25, thick + D * 3.);
		
	vec2 d1 = op_add(groove, dead_wax);
	vec2 d2 = op_add(label, logo);
	vec2 d3 = op_add(d1, d2);
	vec2 d4 = vec2(
		op_sub(d3.x, center_hole),
		d3.y);
		
	return d4;
}

vec3 sdf_normal(_in(vec3) p)
{
	float dt = 0.001;
	vec3 x = vec3(dt, 0, 0);
	vec3 y = vec3(0, dt, 0);
	vec3 z = vec3(0, 0, dt);
	return normalize(vec3(
		sdf(p + x).r - sdf(p - x).r,
		sdf(p + y).r - sdf(p - y).r,
		sdf(p + z).r - sdf(p - z).r
	));
}

vec3 sdf_ao(_in(hit_t) hit)
{
	const float dt = .5;
	const int steps = 5;
	float d = 0.;
	float occlusion = 0.;
	
	for (float i = 1.; i <= float(steps); i += 1.) {
		vec3 p = hit.origin + dt * i * hit.normal;
		d = sdf (p).x;
		
		occlusion += 1. / pow(2., i) * (dt * i - d);
	}
	
	float c = 1. - clamp(occlusion, 0., 1.);
	return vec3 (c, c, c);
}

_constant(vec3) sun_dir = normalize (vec3 (1, 0, 2));

vec3 illuminate(
	_in(vec3) eye,
	_in(hit_t) hit,
	_in(float) ao
) {
#if 0
	return vec3(hit.normal);
#endif

	vec3 V = normalize(eye - hit.origin); // view direction
	vec3 L = sun_dir;
	vec3 accum = vec3(0, 0, 0);
	material_t mat = get_material(hit.material_id);

#if 0
		accum += illum_blinn_phong(V, L, hit, mat);
#else
		accum += illum_cook_torrance(V, L, hit, mat);
#endif

	return accum;
}

vec3 render(
	_in(ray_t) ray,
	_in(vec3) point_cam
){
	const int steps = 70;
	const float end = 20.;

	float t = 0.;
	for (int i = 0; i < steps; i++) {
		vec3 p = ray.origin + ray.direction * t;
		vec2 d = sdf(p);

		if (t > end) break;
		if (d.x < .005) {
			hit_t h = _begin(hit_t)
				t, // ray length at impact
				int(d.y), // material id
				sdf_normal(p),
				p // point of impact				
			_end;

			float ao = 1.;//sdf_ao (h).x;          
			
			return illuminate(ray.origin, h, ao);
		}

		t += d.x;
	}

	return background(ray);
}

#define FOV 1. // 45 degrees
#include "main.h"