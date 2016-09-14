#include "def.h"
#include "util.h"
#include "sdf.h"
#include "material.h"
#include "util_optics.h"
#include "light.h"

#include "noise_iq.h"
#include "fbm.h"
#define noise_func ((noise_iq(p) * 2. - 1.))
DECL_FBM_FUNC(fbm, 4, noise_func)

// ----------------------------------------------------------------------------
// Vinyl disk animation
// ----------------------------------------------------------------------------

vec3 background(_in(ray_t) ray)
{
	return vec3(.6, .6, .6);
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
#if 1
	eye = vec3(0, 6, 7);
	look_at = vec3(0, -3, 0);
#else
	eye = vec3(4, 3, 3);
	look_at = vec3(4, 0, 0);
#endif
}

_mutable(mat3) platter_rot;

float sdf_logo(_in(vec3) pos, _in(float) thick)
{
	vec3 b = vec3(.25, thick, 1.2);
	vec3 d = vec3(.7, 0, 0);
	
	vec3 p = mul(pos, rotate_around_y(30.));
	float v1 = sd_box(p - d, b);
	
	p = mul(pos, rotate_around_y(-30.));
	float v2 = sd_box(p + d, b);
	
	float x = sd_box(pos, vec3(1.5, thick, 1.35));
	float v = op_add(v1, v2);
	return op_intersect(v, x);
}

vec2 sdf(_in(vec3) pos)
{
	const float thick = .1;
	vec3 p = mul(pos, platter_rot);
	
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
		sdf_logo(p, thick - .05),
		mat_logo);
	float center_hole =
		sd_y_cylinder(p, .25, thick * 4.);
		
	vec2 d1 = op_add(groove, dead_wax);
	vec2 d2 = op_add(label, logo);
	vec2 d3 = op_add(d1, d2);
	vec2 d4 = vec2(
		op_sub(d3.x, center_hole),
		d3.y);
#if 0		
	vec2 debug = vec2(
		sd_sphere(p + vec3(4, 0, 0), 1),
		mat_dead_wax);
#endif
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

float saw(_in(float) x)
{
	return x - floor(x);
}

float pulse(_in(float) x)
{
	return saw(x + .5) - saw(x);
}

vec3 illuminate(
	_in(vec3) eye,
	_inout(hit_t) hit
){
#if 0
	return vec3(hit.normal);
#endif

	const vec3 L = normalize(vec3(2, 4, -3));
	const vec3 V = normalize(eye - hit.origin); // view direction

	material_t mat = get_material(hit.material_id);

	if (hit.material_id == mat_groove ||
	hit.material_id == mat_dead_wax) {
		float r = length(hit.origin);
		vec3 B = hit.origin / r;
		vec3 N = vec3(0, 1, 0);
		if (hit.material_id == mat_groove) {
			float rr = r + .075 * fbm(hit.origin * 4.07, 2.08, .5, .5);
			float s = pulse(rr * 6.);
			float ss = pulse(rr * 12.);
			//N.y *= clamp(s, -1., 1.);
			if (s > 0.) {
				N = normalize(N + B);
				if (ss > 0.) {
					N = reflect(N, vec3(0, 1, 0));
				}
			//	N.y = 
			}
		}
		//return N;
		vec3 T = cross(B, N);
		
		vec3 H = normalize(V + L);
		float dotLN = dot(L, N);
		
		const float ro_diff = 1.;
		const float ro_spec = .0625;
		const float a_x = .025;
		const float a_y = .5;
		
		vec3 diffuse = mat.base_color *
			(ro_diff / PI) *
			max(0., dotLN);
			
		float spec_a = ro_spec /
			sqrt(dotLN * dot(V, N));
			
		float spec_b = 1. /
			(4. * PI * a_x * a_y);

		float ht = dot(H, T) / a_x;
		float hb = dot(H, B) / a_y;
		float spec_c = -2. *
			(ht * ht + hb * hb) /
			(1. + dot(H, N));

		vec3 specular = vec3(1, 1, 1) *
			spec_a * spec_b * exp(spec_c);
			
		return diffuse + specular;
	} else {
		hit.normal = sdf_normal(hit.origin);
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
	const int steps = 20;
	const float end = 20.;
	platter_rot = rotate_around_y(u_time * 20.);

	float t = 0.;
	for (int i = 0; i < steps; i++) {
		vec3 p = ray.origin + ray.direction * t;
		vec2 d = sdf(p);

		if (t > end) break;
		if (d.x < .005) {
			hit_t h = _begin(hit_t)
				t, // ray length at impact
				int(d.y), // material id
				vec3(0, 1, 0), // normal
				mul(p, platter_rot) // point of impact				
			_end;    
			
			return illuminate(ray.origin, h);
		}

		t += d.x;
	}

	return background(ray);
}

#define FOV 1. // 45 degrees
#include "main.h"