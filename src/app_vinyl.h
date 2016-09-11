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
	setup_mat(materials[mat_debug],
		vec3(1, 1, 1), .0, .0);
	setup_mat(materials[mat_groove],
		vec3(.01, .01, .01), .0, .013);
	setup_mat(materials[mat_dead_wax],
		vec3(.05, .05, .05), .0, .005);
	setup_mat(materials[mat_label],
		vec3(.5, .5, .0), .0, .5);
	setup_mat(materials[mat_logo],
		vec3(0, 0, .7), .0, .5);
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
	
	vec3 p = mul(pos, rotate_around_x(90.));
	p = mul(p, rotate_around_y(u_time * 200.));
	const float thick = .1;
	
	vec2 groove = vec2(
		sd_y_cylinder(p, 6., thick),
		mat_groove);
	vec2 dead_wax = vec2(
		sd_y_cylinder(p, 2.5, thick),
		mat_dead_wax);
	vec2 label = vec2(
		sd_y_cylinder(p, 2., thick),
		mat_label);
	vec2 logo = vec2(
		op_add(
			sd_box(p - vec3(0, 0, .9), vec3(1., thick, .5)),
			sd_box(p + vec3(0, 0, .8), vec3(1., thick, .25))),
		mat_logo);
	float center_hole =
		sd_y_cylinder(p, .25, thick * 4.);
		
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
#if 0
	return vec3(0,0,1);
#endif
	
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

vec3 illuminate(
	_in(vec3) eye,
	_in(hit_t) hit
){
#if 0
	return vec3(hit.normal);
#endif

	const vec3 L = normalize(vec3(1, 0, 2));
	const vec3 V = normalize(eye - hit.origin); // view direction

	material_t mat = get_material(hit.material_id);

	if (hit.material_id == mat_groove ||
	hit.material_id == mat_dead_wax) {
		vec3 B = normalize(hit.origin);
		vec3 T = cross(B, hit.normal);
		
		vec3 H = normalize(V + L);
		float dotLN = dot(L, hit.normal);
		
		const float ro_diff = .1;
		const float ro_spec = .33;
		const float a_x = .05;
		const float a_y = .16;
		
		vec3 diffuse = mat.base_color *
			(ro_diff / PI) *
			max(0., dotLN);
			
		float spec_a = ro_spec /
			sqrt(dotLN * dot(V, hit.normal));
			
		float spec_b = 1. /
			(4. * PI * a_x * a_y);

		float ht = dot(H, T) / a_x;
		float hb = dot(H, B) / a_y;
		float spec_c = -2. *
			(ht * ht + hb * hb) /
			(1. + dot(H, hit.normal));

		vec3 specular = vec3(1, 1, 1) *
			spec_a * spec_b * exp(spec_c);
			
		return diffuse + specular;
	} else {
#if 1
		return illum_blinn_phong(V, L, hit, mat);
#else
		return illum_cook_torrance(V, L, hit, mat);
#endif
	}
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
			
			return illuminate(ray.origin, h);
		}

		t += d.x;
	}

	return background(ray);
}

#define FOV 1. // 45 degrees
#include "main.h"