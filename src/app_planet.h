#include "def.h"
#include "util.h"
#include "intersect.h"

#include "noise_iq.h"
#define noise(x) noise_iq(x)
#include "fbm.h"

// ----------------------------------------------------------------------------
// Planet
// ----------------------------------------------------------------------------

_constant(sphere_t) planet = _begin(sphere_t)
vec3(0, 0, 0), 1., 0
_end;

DECL_FBM_FUNC(fbm_terr, 4, .5)
DECL_FBM_FUNC(fbm_terr_nrm, 7, .5)
DECL_FBM_FUNC(fbm_cloud, 5, .5)

#define max_height .35
#define max_ray_dist (max_height * 4.)

#define TERR_STEPS 120
#define TERR_EPS .005
#define TERR_FLAT_BIAS .35
#define TERR_FBM_FREQ 1.9870725
#define TERR_FBM_LAC 2.023674

#define CLOUDS
#define CLD_STEPS 75
#define t_cld_step (max_ray_dist / float(CLD_STEPS))
#define CLD_COVERAGE .575675 // higher=less clouds
#define CLD_FUZZY .0335 // higher=fuzzy, lower=blockier
#define CLD_ABSORBTION 33.93434
#define CLD_TOP .9
#define CLD_LOW .5
#define CLD_MED (CLD_LOW + (CLD_TOP - CLD_LOW) / 2.)

vec3 background(
	_in(ray_t) eye
){
	_constant(vec3) sun_color = vec3(1., .9, .55);
	float sun_amount = dot(eye.direction, vec3(0, 0, 1));

	vec3 sky = mix(
		vec3(.0, .05, .2),
		vec3(.15, .3, .4),
		1.0 - eye.direction.y);
	sky += sun_color * min(pow(sun_amount, 30.0) * 5.0, 1.0);
	sky += sun_color * min(pow(sun_amount, 10.0) * .6, 1.0);

	return sky;
}

void setup_scene()
{
}

void setup_camera(
	_inout(vec3) eye,
	_inout(vec3) look_at
){
	eye = vec3(0, 0, -2.5);
	look_at = vec3(0, 0, 2);
}

// ----------------------------------------------------------------------------
// Clouds
// ----------------------------------------------------------------------------

struct cloud_drop_t {
	vec3 origin;
	vec3 pos;
	float T;
	vec3 C;
	float alpha;
};

_mutable(cloud_drop_t) cloud;

cloud_drop_t start_cloud(
	_in(vec3) origin
){
	cloud_drop_t c = _begin(cloud_drop_t)
		origin, origin, 1., vec3 (0., 0., 0.), 0.
	_end;
	return c;
}

void clouds_map(
	_inout(cloud_drop_t) cloud,
	_in(float) height,
	_in(float) t_step
){
	float dens = fbm_cloud(
		cloud.pos * 2.2343 + vec3(1.35, 3.35, 2.67),
		2.9760);

	dens *= smoothstep(CLD_COVERAGE, CLD_COVERAGE + CLD_FUZZY, dens);

	dens *= band(CLD_LOW, CLD_MED, CLD_TOP, height);

	float T_i = exp(-CLD_ABSORBTION * dens * t_step);
	cloud.T *= T_i;
	cloud.C +=
		cloud.T * (exp(height) / .055)
		* dens * t_step;
	cloud.alpha += (1. - T_i) * (1. - cloud.alpha);
}

void clouds_march(
	_in(ray_t) eye,
	_inout(cloud_drop_t) cloud,
	_in(float) max_travel,
	_in(mat3) rot
){
	float t = 0.;

	for (int i = 0; i < CLD_STEPS; i++) {
		if (t > max_travel || cloud.alpha >= 1.) return;

		vec3 o = cloud.origin + t * eye.direction;
		cloud.pos = mul(rot, o - planet.origin);

		float h = (length(cloud.pos) - planet.radius) / max_height;
		t += t_cld_step;
		clouds_map(cloud, h, t_cld_step);
	}
}

void clouds_shadow_march(
	_in(vec3) dir,
	_inout(cloud_drop_t) cloud,
	_in(mat3) rot
){
	const int steps = 5;
	const float t_step = max_height / float(steps);
	float t = 0.;

	for (int i = 0; i < steps; i++) {
		vec3 o = cloud.origin + t * dir;
		cloud.pos = mul(rot, o - planet.origin);

		float h = (length(cloud.pos) - planet.radius) / max_height;
		t += t_step;
		clouds_map(cloud, h, t_step);
	}
}

// ----------------------------------------------------------------------------
// Terrain
// ----------------------------------------------------------------------------

vec2 sdf_terrain_map(_in(vec3) pos)
{
	float h = fbm_terr(pos * TERR_FBM_FREQ, TERR_FBM_LAC);
	float n = smoothstep(TERR_FLAT_BIAS, 1., h);
	return vec2(length(pos) - planet.radius - n * max_height, n / max_height);
}

vec2 sdf_terrain_map_detail(_in(vec3) pos)
{
	float h = fbm_terr_nrm(pos * TERR_FBM_FREQ, TERR_FBM_LAC);
	float n = smoothstep(TERR_FLAT_BIAS, 1., h);
	return vec2(length(pos) - planet.radius - n * max_height, n);
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

	// fill light 2 - ambient (reversed key with no height)
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
//TODO: moves as much as possible in #define's
//TODO: paste final material values
		// current terrain height at position
	float h = df.y;

	vec3 normal = sdf_terrain_normal(pos);
	vec3 w_normal = normalize(pos);
	//return abs(normal);
	float N = abs(normal.y);
	//dot(normal, normalize(pos));

// materials
	vec3 c_water = srgb_to_linear(vec3(38, 94, 179) / 256.);
	vec3 c_grass = srgb_to_linear(vec3(84, 102, 42) / 256.);
	vec3 c_beach = srgb_to_linear(vec3(109, 115, 98) / 256.);
	vec3 c_rock1 = vec3(0.08, 0.05, 0.03);
	vec3 c_rock2 = vec3(0.10, 0.09, 0.08);
	vec3 c_snow = vec3(1., 1., 1.) * .4;

	// limits
	const float l_water = .05;
	const float l_shore = .1;
	const float l_rock = .211;

	vec3 rock_strata = mix(
		c_rock1, c_rock2,
		smoothstep(l_rock, 1.,
			cos(h * 45.47 * fbm(pos * 17.24, 4.37))));
	float s = smoothstep(.5, 1., h);
	vec3 rock = mix(
		rock_strata, c_snow,
		//smoothstep(l_rock, 1., h));
		smoothstep(1. - .4*s, 1. - .1*s, N));

	vec3 shoreline = mix(
		c_beach, rock,
		smoothstep(l_shore, l_rock, h));

	// lights
	vec3 L = mul(local_xform, vec3(1, 0, 0));
		//mul(local_xform, mul(rotate_around_y(sin(u_time) * 50.0), normalize(vec3(-1, 1, 0))));
	shoreline *= setup_lights(L, normal);
	c_water *= setup_lights(L, w_normal);

	return mix(
		c_water, shoreline,
		smoothstep(l_water, l_shore, h));
}

vec3 render(
	_in(ray_t) eye,
	_in(vec3) point_cam
){
	mat3 rot = rotate_around_x(u_time * 8.);
	mat3 rot_cloud = rotate_around_x(u_time * -8.);

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

		t += df.x *.67;
	}

#ifdef CLOUDS
	cloud = start_cloud(hit.origin);
	clouds_march(eye, cloud, max_cld_ray_dist, rot_cloud);
#endif
	
	if (df.x < TERR_EPS) {
		vec3 c_terr = illuminate(pos, eye.direction, rot, df);
		vec3 c_cld = cloud.C;
		float alpha = cloud.alpha;
		float shadow = 1.;

#ifdef CLOUDS // clouds ground shadows
		pos = mul(transpose(rot), pos);
		cloud = start_cloud(pos);
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