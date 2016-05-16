#include "def.h"
#include "util.h"
#include "intersect.h"
#include "volumetric.h"

#include "noise_iq.h"
#define noise(x) noise_iq(x)

#include "fbm.h"

// ----------------------------------------------------------------------------
// Planet
// ----------------------------------------------------------------------------
_constant(sphere_t) planet = _begin(sphere_t)
	vec3(0, 0, 0), 1., 0
_end;

#define max_height .4
#define max_ray_dist (max_height * 4.)

vec3 background(
	_in(ray_t) eye
){
#if 1
	return vec3(.15, .3, .4);
#else
	_constant(vec3) sun_color = vec3(1., .9, .55);
	float sun_amount = dot(eye.direction, vec3(0, 0, 1));

	vec3 sky = mix(
		vec3(.0, .05, .2),
		vec3(.15, .3, .4),
		1.0 - eye.direction.y);
	sky += sun_color * min(pow(sun_amount, 30.0) * 5.0, 1.0);
	sky += sun_color * min(pow(sun_amount, 10.0) * .6, 1.0);

	return sky;
#endif
}

void setup_scene()
{
}

void setup_camera(
	_inout(vec3) eye,
	_inout(vec3) look_at
){
#if 0
	eye = vec3(.0, 0, -1.93);
	look_at = vec3(-.1, .9, 2);
#else
	eye = vec3(0, 0, -2.5);
	look_at = vec3(0, 0, 2);
#endif
}

// ----------------------------------------------------------------------------
// Clouds
// ----------------------------------------------------------------------------
#define CLOUDS

#define anoise (abs(noise(p) * 2. - 1.))
DECL_FBM_FUNC(fbm_clouds, 4, anoise)

#define vol_coeff_absorb 30.034
_mutable(volume_sampler_t) cloud;

float illuminate_volume(
	_inout(volume_sampler_t) cloud,
	_in(vec3) V,
	_in(vec3) L
){
	return exp(cloud.height) / .055;
}

void clouds_map(
	_inout(volume_sampler_t) cloud,
	_in(float) t_step
){
	float dens = fbm_clouds(
		cloud.pos * 3.2343 + vec3(.35, 13.35, 2.67),
		2.0276, .5, .5);

	#define cld_coverage .29475675 // higher=less clouds
	#define cld_fuzzy .0335 // higher=fuzzy, lower=blockier
	dens *= smoothstep(cld_coverage, cld_coverage + cld_fuzzy, dens);

	dens *= band(.2, .35, .65, cloud.height);

	integrate_volume(cloud,
		cloud.pos, cloud.pos, // unused dummies 
		dens, t_step);
}

void clouds_march(
	_in(ray_t) eye,
	_inout(volume_sampler_t) cloud,
	_in(float) max_travel,
	_in(mat3) rot
){
	const int steps = 75;
	const float t_step = max_ray_dist / float(steps);
	float t = 0.;

	for (int i = 0; i < steps; i++) {
		if (t > max_travel || cloud.alpha >= 1.) return;
			
		vec3 o = cloud.origin + t * eye.direction;
		cloud.pos = mul(rot, o - planet.origin);

		cloud.height = (length(cloud.pos) - planet.radius) / max_height;
		t += t_step;
		clouds_map(cloud, t_step);
	}
}

void clouds_shadow_march(
	_in(vec3) dir,
	_inout(volume_sampler_t) cloud,
	_in(mat3) rot
){
	const int steps = 5;
	const float t_step = max_height / float(steps);
	float t = 0.;

	for (int i = 0; i < steps; i++) {
		vec3 o = cloud.origin + t * dir;
		cloud.pos = mul(rot, o - planet.origin);

		cloud.height = (length(cloud.pos) - planet.radius) / max_height;
		t += t_step;
		clouds_map(cloud, t_step);
	}
}

// ----------------------------------------------------------------------------
// Terrain
// ----------------------------------------------------------------------------
#define TERR_STEPS 120
#define TERR_EPS .005
#define rnoise (1. - abs(noise(p) * 2. - 1.))

DECL_FBM_FUNC(fbm_terr, 3, noise(p))
DECL_FBM_FUNC(fbm_terr_r, 3, rnoise)

DECL_FBM_FUNC(fbm_terr_normals, 7, noise(p))
DECL_FBM_FUNC(fbm_terr_r_normals, 7, rnoise)

vec2 sdf_terrain_map(_in(vec3) pos)
{
	float h0 = fbm_terr(pos * 2.0987, 2.0244, .454, .454);
	float n0 = smoothstep(.35, 1., h0);

	float h1 = fbm_terr_r(pos * 1.50987 + vec3(1.9489, 2.435, .5483), 2.0244, .454, .454);
	float n1 = smoothstep(.6, 1., h1);
	
	float n = n0 + n1;
	
	return vec2(length(pos) - planet.radius - n * max_height, n / max_height);
}

vec2 sdf_terrain_map_detail(_in(vec3) pos)
{
	float h0 = fbm_terr_normals(pos * 2.0987, 2.0244, .454, .454);
	float n0 = smoothstep(.35, 1., h0);

	float h1 = fbm_terr_r_normals(pos * 1.50987 + vec3(1.9489, 2.435, .5483), 2.0244, .454, .454);
	float n1 = smoothstep(.6, 1., h1);

	float n = n0 + n1;

	return vec2(length(pos) - planet.radius - n * max_height, n / max_height);
}

vec3 sdf_terrain_normal(_in(vec3) p)
{
#define F(t) sdf_terrain_map_detail(t).x
	vec3 dt = vec3(0.001, 0, 0);

	return normalize(vec3(
		F(p + dt.xzz) - F(p - dt.xzz),
		F(p + dt.zxz) - F(p - dt.zxz),
		F(p + dt.zzx) - F(p - dt.zzx)
	));
#undef F
}

// ----------------------------------------------------------------------------
// Lighting
// ----------------------------------------------------------------------------
vec3 setup_lights(
	_in(vec3) L,
	_in(vec3) normal
){
	vec3 diffuse = vec3(0, 0, 0);

	// key light
	vec3 c_L = vec3(7, 5, 3);
	diffuse += max(0., dot(L, normal)) * c_L;

	// fill light 1 - faked hemisphere
	float hemi = clamp(.25 + .5 * normal.y, .0, 1.);
	diffuse += hemi * vec3(.4, .6, .8) * .2;

	// fill light 2 - ambient (reversed key)
	float amb = clamp(.12 + .8 * max(0., dot(-L, normal)), 0., 1.);
	diffuse += amb * vec3(.4, .5, .6);

	return diffuse;
}

vec3 illuminate(
	_in(vec3) pos,
	_in(vec3) eye,
	_in(mat3) local_xform,
	_in(vec2) df
){
	// current terrain height at position
	float h = df.y;
	//return vec3 (h);

	vec3 w_normal = normalize(pos);
#define LIGHT
#ifdef LIGHT
	vec3 normal = sdf_terrain_normal(pos);
	float N = dot(normal, w_normal);
#else
	float N = w_normal.y;
#endif

	// materials
	#define c_water vec3(.015, .110, .455)
	#define c_grass vec3(.086, .132, .018)
	#define c_beach vec3(.153, .172, .121)
	#define c_rock  vec3(.080, .050, .030)
	#define c_snow  vec3(.600, .600, .600)

	// limits
	#define l_water .05
	#define l_shore .17
	#define l_grass .211
	#define l_rock .351

	float s = smoothstep(.4, 1., h);
	vec3 rock = mix(
		c_rock, c_snow,
		smoothstep(1. - .3*s, 1. - .2*s, N));

	vec3 grass = mix(
		c_grass, rock,
		smoothstep(l_grass, l_rock, h));
		
	vec3 shoreline = mix(
		c_beach, grass,
		smoothstep(l_shore, l_grass, h));

	vec3 water = mix(
		c_water / 2., c_water,
		smoothstep(0., l_water, h));

#ifdef LIGHT
	vec3 L = mul(local_xform, normalize(vec3(1, 1, 0)));
	shoreline *= setup_lights(L, normal);
	vec3 ocean = setup_lights(L, w_normal) * water;
#else
	vec3 ocean = water;
#endif
	
	return mix(
		ocean, shoreline,
		smoothstep(l_water, l_shore, h));
}

// ----------------------------------------------------------------------------
// Rendering
// ----------------------------------------------------------------------------
vec3 render(
	_in(ray_t) eye,
	_in(vec3) point_cam
){
	mat3 rot_y = rotate_around_y(27.);
	mat3 rot = mul(rotate_around_x(u_time * -12.), rot_y);
	mat3 rot_cloud = mul(rotate_around_x(u_time * 8.), rot_y);

	sphere_t atmosphere = planet;
	atmosphere.radius += max_height;

	hit_t hit = no_hit;
	intersect_sphere(eye, atmosphere, hit);
	if (hit.material_id < 0) {
		return background(eye);
	}

	float t = 0.;
	vec2 df = vec2(1, max_height);
	vec3 pos;
	float max_cld_ray_dist = max_ray_dist;
	
	for (int i = 0; i < TERR_STEPS; i++) {
		if (t > max_ray_dist) break;
		
		vec3 o = hit.origin + t * eye.direction;
		pos = mul(rot, o - planet.origin);

		df = sdf_terrain_map(pos);

		if (df.x < TERR_EPS) {
			max_cld_ray_dist = t;
			break;
		}

		t += df.x * .4567;
	}

#ifdef CLOUDS
	cloud = begin_volume(hit.origin, vol_coeff_absorb);
	clouds_march(eye, cloud, max_cld_ray_dist, rot_cloud);
#endif
	
	if (df.x < TERR_EPS) {
		vec3 c_terr = illuminate(pos, eye.direction, rot, df);
		vec3 c_cld = cloud.C;
		float alpha = cloud.alpha;
		float shadow = 1.;

#ifdef CLOUDS // clouds ground shadows
		pos = mul(transpose(rot), pos);
		cloud = begin_volume(pos, vol_coeff_absorb);
		vec3 local_up = normalize(pos);
		clouds_shadow_march(local_up, cloud, rot_cloud);
		shadow = mix(.7, 1., step(cloud.alpha, 0.33));
#endif

		return mix(c_terr * shadow, c_cld, alpha);
	} else {
		return mix(background(eye), cloud.C, cloud.alpha);
	}
}

#define FOV tan(radians(30.))
#include "main.h"