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

#define TERR_STEPS 120
#define TERR_EPS .005
#define max_height .35
#define max_ray_dist (max_height * 4.)

#define CLD_STEPS 75
#define t_cld_step (max_ray_dist / float(CLD_STEPS))

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

struct cloud_drop_t {
	float t;
	vec3 origin;
	vec3 pos;
	float h;
	float T;
	vec3 C;
	float alpha;
};

cloud_drop_t start_cloud(
	_in(vec3) origin
){
	cloud_drop_t c = _begin(cloud_drop_t)
		0., origin, origin, 0., 1., vec3 (0., 0., 0.), 0.
	_end;
	return c;
}

_mutable(cloud_drop_t) cloud;

DECL_FBM_FUNC(fbm_cloud, 5, .5)

void clouds_map(
	_inout(cloud_drop_t) cloud,
	_in(float) t_step
){
	float H = exp(cloud.h);
	float dens = fbm_cloud(cloud.pos * 1.2343
		+ vec3(1.35, 3.35, 2.67), 2.9760);

	const float coverage = .4575675; // higher=less clouds
	const float fuzziness = .0335; // higher=fuzzy, lower=blockier
	//dens *= step(coverage, dens);
	dens *= smoothstep(coverage, coverage + fuzziness, dens);

	// these 2 are identical ways to "band" vertically
	//TODO: understand why the exp thing really works
	//dens *= band(.2, .4, .6, exp(h) / 4.);
	dens *= 1. - smoothstep(.4, .5, H / 4.);

	const float absorbtion = 33.93434;
	float T_i = exp(-absorbtion * dens * t_step);
	cloud.T *= T_i;
	cloud.C +=
		cloud.T * (H / .055)
		* dens * t_step;
	cloud.alpha += (1. - T_i) * (1. - cloud.alpha);
}

void clouds_march(
	_in(ray_t) eye,
	_inout(cloud_drop_t) cloud,
	_in(float) max_travel,
	_in(mat3) rot
){
	for (int i = 0; i < CLD_STEPS; i++) {
		if (cloud.t > max_travel || cloud.alpha >= 1.) return;

		vec3 o = cloud.origin + cloud.t * eye.direction;
		cloud.pos = mul(rot, o - planet.origin);

		cloud.h = (length(cloud.pos) - planet.radius) / max_height;
		cloud.t += t_cld_step;

		clouds_map(cloud, t_cld_step);
	}
}

DECL_FBM_FUNC(fbm_terr, 4, .5)
DECL_FBM_FUNC(fbm_terr_nrm, 7, .5)

#define TERR_FLAT_BIAS .35
#define TERR_FBM_FREQ 1.9870725
#define TERR_FBM_LAC 2.023674

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

vec3 render_planet(
	_in(ray_t) eye
){
	mat3 rot = rotate_around_x(u_time * 16.);
	mat3 rot2 = rotate_around_x(u_time * -32.);

	sphere_t atmosphere = planet;
	atmosphere.radius += max_height;

	hit_t hit = no_hit;
	intersect_sphere(eye, atmosphere, hit);
	if (hit.material_id < 0) {
		return background(eye);
	}

#if 0 // test with checkboard pattern
	vec3 d = mul(rot, hit.normal);
#ifdef HLSL
#define atan(y, x) atan2(x, y)
#endif
	float u = .5 + atan(d.z, d.x) / (2. * PI);
	float v = .5 - asin(d.y) / (1. * PI);
	float n = checkboard_pattern(vec2(u, v), 20.);
	return vec3(n, n, n);
#endif

	float t = 0.;
	vec2 df = vec2(1, max_height);
	vec3 c_out = vec3(0, 0, 0);

	cloud = start_cloud (hit.origin);
	float max_cld_ray_dist = max_ray_dist;
	
	for (int i = 0; i < TERR_STEPS; i++) {
		if (t > max_ray_dist) break;
		
		vec3 o = hit.origin + t * eye.direction;
		vec3 p = mul(rot, o - planet.origin);

		df = sdf_terrain_map(p);

		if (df.x < TERR_EPS) {
			float h = df.y;
			
			vec3 normal = sdf_terrain_normal(p);
			//return abs(normal);
			float N = abs(normal.y);
				//dot(normal, normalize(p));

			// colours TODO: paste in final values
			vec3 c_water = srgb_to_linear(vec3( 38,  94, 179) / 256.);
			vec3 c_grass = srgb_to_linear(vec3( 84, 102,  42) / 256.);
			vec3 c_beach = srgb_to_linear(vec3(109, 115,  98) / 256.);
			vec3 c_rock1 = srgb_to_linear(vec3(103, 109, 111) / 256.);
			vec3 c_rock2 = srgb_to_linear(vec3(100, 107,  98) / 256.);
			vec3 c_snow  = vec3(1, 1, 1);

			// limits TODO: vary with fbm
			const float l_water = .05;
			const float l_shore = .1;
			const float l_rock = .211;
			
			// light
			const vec3 L = vec3(1, 0, 0);
			const vec3 c_L = vec3(1, 1, 1) * 3.;
			const vec3 c_amb = vec3(.1, .1, .1);

			//vec3 rock = mix(
			//	c_rock1, c_rock2,
			//	smoothstep(l_rock, 1.,
			//		cos(hs * 85.47 * fbm(p * 17.24, 4.37))));
			float s = smoothstep(.5, 1., h);
			vec3 rock = mix(
				c_rock2, c_snow,
				smoothstep(1. - .4*s, 1. - .1*s, N));

			vec3 shoreline = mix(
				c_beach, rock,
				smoothstep(l_shore, l_rock, h));
			shoreline *= max(0., dot(L, normal)) * c_L;
			
			c_out = mix(
				c_water, shoreline,
				smoothstep (l_water, l_shore, h));
			
			max_cld_ray_dist = t;
			break;
		}

		t += df.x *.67;
	}

	clouds_march(eye, cloud, max_cld_ray_dist, rot2);
	
	if (df.x < TERR_EPS) {
		return mix(c_out, cloud.C, cloud.alpha);
	} else {
		return mix(background(eye), cloud.C, cloud.alpha);
	}
}

vec3 render(
	_in(ray_t) eye,
	_in(vec3) point_cam
){
	return render_planet(eye);
}

#define FOV tan(radians(30.))
#include "main.h"