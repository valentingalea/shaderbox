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
	return vec3(1, 1, 1);
}

#define mat_groove 1
#define mat_dead_wax 2
#define mat_label 3
#define mat_logo 4
#define mat_shiny 5

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
	setup_mat(materials[mat_shiny],
		vec3(.7, .7, .7), 1., .01);
}

void setup_camera(_inout(vec3) eye, _inout(vec3) look_at)
{
#if 1
	eye = vec3(0, 5.75, 6.75);
	look_at = vec3(0, -2.5, 0);
#else
	eye = vec3(0, 0, 15);
	look_at = vec3(0, 0, 0);
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

vec2 sdf_platter(_in(vec3) pos)
{
	const float thick = .1;
	vec3 p = mul(pos, platter_rot);

	vec2 lead_in = vec2(
		sd_y_cylinder(p, 6., thick - .05),
		mat_dead_wax);
	vec2 groove = vec2(
		sd_y_cylinder(p, 5.9, thick),
		mat_groove);
	vec2 dead_wax = vec2(
		sd_y_cylinder(p, 3., thick),
		mat_dead_wax);
	vec2 label = vec2(sd_y_cylinder(p, 2., thick),
		mat_label);
	vec2 logo = vec2(
		sdf_logo(p, thick - .0175),
		mat_logo);
	//float center_hole =
	//	sd_y_cylinder(p, .25, thick * 4.);
	float spc = sd_y_cylinder(p, .10, .6);
	float sps = sd_sphere(p - vec3(0, .3, 0), .10);
	vec2 spindle = vec2(
		op_add(spc, sps), mat_shiny);
	
	vec2 d0 = op_add(groove, lead_in);
	vec2 d1 = op_add(d0, dead_wax);
	vec2 d2 = op_add(label, logo);
	vec2 d3 = op_add(d1, d2);
	vec2 d4 = op_add(d3, spindle);
	//	op_sub(d3.x, center_hole),
	//	d3.y);

	float defect1 = sd_sphere(p + vec3(6.05, 0, 0), .1);
	float defect2 = sd_sphere(p + vec3(-6.05, 0, 0), .1);
	float defect = op_add(defect1, defect2);

	return vec2(op_sub(d4.x, defect), d4.y);
}

// Distance to line segment between <a> and <b>, used for fCapsule() version 2below
float fLineSegment(vec3 p, vec3 a, vec3 b) {
	vec3 ab = b - a;
	float t = clamp(dot(p - a, ab) / dot(ab, ab), 0., 1.);
	return length((ab*t + a) - p);
}

// Capsule version 2: between two end points <a> and <b> with radius r 
float sd_capsule(vec3 p, vec3 a, vec3 b, float r) {
	return fLineSegment(p, a, b) - r;
}

vec2 sdf(_in(vec3) p)
{
	vec2 plat = sdf_platter(p);

	vec3 base_p = vec3(-7, 0, -5);
	float platter = sd_y_cylinder(p, 6.25, 1.);
	float base_0 = sd_y_cylinder(p - base_p, 3., .25);
	float base_1 = op_sub(base_0, platter);
	float base_2 = sd_y_cylinder(p - base_p, 1.15, 1.05);
	float base_12 = op_add(base_1, base_2);
	vec2 base_a = vec2(base_12, mat_shiny);
	vec2 base_b = vec2(sd_y_cylinder(p - base_p, 0.5, 2.5), mat_shiny);
	vec2 base = op_add(base_a, base_b);

	const float D = .15;
	const float H = 1.;
	vec3 a1 = vec3(-6, H, -3);
	vec3 a11 = vec3(-4.25, H, 2);
	vec3 a2 = vec3(-4.1, H, 2.45);
	vec3 a33 = vec3(-3.5, H, 3);
	vec3 a3 = vec3(-2, H, 4);
	float arm1 = sd_capsule(p, base_p + vec3(-1, H, -2), a1, D);
	float arm2 = sd_capsule(p, a1, a11, D);
	float arm3 = sd_capsule(p, a33, a3, D);
	vec2 armb = sd_bezier(a11, a2, a33, p, D);
	vec2 arm = vec2(
		op_add(op_add(op_add(arm1, arm2), arm3), armb.x),
		mat_shiny);

	vec3 ang = vec3(20, 35, 45);
#ifdef HLSL
	ang *= -1.;
#endif
	vec3 pp = mul(p - a3,
		mul(mul(rotate_around_y(ang.y),
		rotate_around_z(ang.z)),
		rotate_around_x(ang.x)));
	vec2 cartrige = vec2(
		sd_box(pp - vec3(.45, 0, 0), vec3(.5, .15, .15)),
		mat_shiny);

	vec2 tonearm = op_add(op_add(base, arm), cartrige);

	return op_add(plat, tonearm);
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

_constant(vec3) sun_dir = //normalize(vec3(2, 4, -3));
	normalize(vec3(-1, 4, -3));

vec3 illuminate(
	_in(vec3) eye,
	_inout(hit_t) hit
){
#if 0
	return vec3(hit.normal);
#endif

	vec3 L = sun_dir;
	vec3 V = normalize(eye - hit.origin); // view direction

	material_t mat = get_material(hit.material_id);

	if (hit.material_id == mat_groove ||
	hit.material_id == mat_dead_wax) {
		hit.origin = mul(hit.origin, platter_rot);
		L = mul(L, platter_rot);
		V = mul(V, platter_rot);
		
		float r = length(hit.origin);
		vec3 B = hit.origin / r;
		vec3 N = vec3(0, 1, 0);
		if (hit.material_id == mat_groove) {
			float rr = r + .07575 *
				noise_iq(hit.origin * 2.456);
				//fbm(hit.origin * 4.07, 2.08, .5, .5);
			float s = pulse(rr * 24.);
			if (s > 0.) {
				N = normalize(N + B);
				N = reflect(N, vec3(0, 1, 0));
			}
		}
		if (hit.material_id == mat_dead_wax) {
			float s = saw(r * 4.);
			N = normalize(N + B * float(s > .9));
		}
		//return N;
		vec3 T = cross(B, N);
		
		vec3 H = normalize(V + L);
		float dotLN = dot(L, N);
		
		const float ro_diff = 1.;
		const float ro_spec = .0725;
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
#if 0
		if (hit.material_id == mat_label || hit.material_id == mat_logo) {
			float r = length(hit.origin);
			vec3 B = hit.origin / r;
			float s = saw(r * .9);
			hit.normal = normalize(hit.normal + B * float(s > .975));
		}
#endif

#if 1
		return illum_blinn_phong(V, L, hit, mat);
#else
		return illum_cook_torrance(V, L, hit, mat);
#endif
	}
}

float sdf_shadow(_in(ray_t) ray)
{
	const int steps = 20;
	const float end = 5.;
	const float penumbra_factor = 16.;
	const float darkest = .05;
	float t = 0.;
	float umbra = 1.;

	for (int i = 0; i < steps; i++) {
		vec3 p = ray.origin + ray.direction * t;
		vec2 d = sdf(p);

		if (t > end) break;
		if (d.x < .005) {
			return darkest;
		}

		t += d.x;
		// from http://iquilezles.org/www/articles/rmshadows/rmshadows.htm
		umbra = min(umbra, penumbra_factor * d.x / t);
	}

	return umbra;
}

vec3 render(
	_in(ray_t) ray,
	_in(vec3) point_cam
){
	const int steps = 80;
	const float end = 40.;

	float rot = u_time * 200.;
#ifdef SHADERTOY
	if (u_mouse.z > 0.) {
		rot = u_mouse.y;
	}
#endif
	platter_rot = mul(
		rotate_around_y(rot), 
		rotate_around_x(sin(u_time) * .1));

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
				p // point of impact				
			_end;    
		
			float sh = 1.;
#if 1
			ray_t sh_ray = _begin(ray_t)
				p + sun_dir * 0.05, sun_dir
				_end;
			sh = sdf_shadow(sh_ray);
#endif   
			return illuminate(ray.origin, h) * sh;
		}

		t += d.x;
	}

	return background(ray);
}

#define FOV 1. // 45 degrees
#include "main.h"