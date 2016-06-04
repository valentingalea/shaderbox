#include "def.h"

#include "../lib/ashima-noise/src/common.glsl"
#include "../lib/ashima-noise/src/classicnoise3d.glsl"
#include "../lib/ashima-noise/src/noise3d.glsl"
#include "../lib/ashima-noise/src/cellular3d.glsl"
#include "noise_worley.h"

#include "fbm.h"
DECL_FBM_FUNC(fbm_perlin, 4, abs(cnoise(p)))
DECL_FBM_FUNC(fbm_simplex, 4, abs(snoise(p)))
DECL_FBM_FUNC(fbm_worley, 3, (1. - cellular(p).r))

DECL_FBM_FUNC_TILE(fbm_worley_tile, 4, (1. - (noise_w(p, L).r + .25)))
DECL_FBM_FUNC_TILE(fbm_perlin_tile, 4, abs(pcnoise(p, L)))

float worley_tex_left(_in(vec3) pos)
{
	float w1 = (1. - (noise_w(pos, 4.).r + .015));
	float w2 = (1. - (noise_w(pos, 8.).r + .015));
	float w3 = (1. - (noise_w(pos, 16.).r + .015));
	return w1 * .625 + w2 * .25 + w3 * .125;
}

float worley_tex_middle(_in(vec3) pos)
{
	float w1 = (1. - (noise_w(pos, 8.).r + .015));
	float w2 = (1. - (noise_w(pos, 16.).r + .015));
	float w3 = (1. - (noise_w(pos, 32.).r + .015));
	return w1 * .625 + w2 * .25 + w3 * .125;
}

float worley_tex_right(_in(vec3) pos)
{
	float w1 = (1. - (noise_w(pos, 24.).r + .015));
	float w2 = (1. - (noise_w(pos, 32.).r + .015));
	float w3 = (1. - (noise_w(pos, 64.).r + .015));
	return w1 * .625 + w2 * .25 + w3 * .125;
}

float worley_fbm(_in(vec3) pos)
{
	float w1 = worley_tex_left(pos);
	float w2 = worley_tex_middle(pos);
	float w3 = worley_tex_right(pos);
	return w1 * .625 + w2 * .25 + w3 * .125;
}

#define SCALE 1.
#define D (.0125 * SCALE)

vec3 plot(
	_in(float) f,
	_in(float) x,
	_in(vec3) color
){
	float y = smoothstep (f-D, f+D, x);
	y *= 1.- y;
	
	return y * color * 5.;
}

void mainImage(
	_out(vec4) fragColor,
#ifdef SHADERTOY
	vec2 fragCoord
#else
	_in(vec2) fragCoord
#endif
){
// go from [0..resolution] to [0..1]
	vec2 t = (fragCoord.xy + .5) / u_res.xy;
#ifdef HLSL
	t.y = 1. - t.y;
#endif

	vec3 col = vec3(0, 0, 0);

#if 1 // 2D
	vec3 pos = vec3(t, 0);
	float n =
		//fbm_worley_tile(pos, 2, 1., .5);
		worley_fbm(pos);

	col += vec3(n, n, n);
#else // 1D

// optional: center around origin by going to [-1, +1] and scale
	t = (t * 2. - 1.) * SCALE;
	
// plot the axes
	col += plot(0., t.y, vec3(1, 1, 1));
	col += plot(t.x, 0., vec3(1, 1, 1));

// optional: animate
	t.x += u_time * SCALE;

// plot custom functions	
	//col += plot(sin(t.x), t.y, vec3(1, 0, 0));
	//col += plot(cos(t.x), t.y, vec3(0, 1, 0));
	//col += plot(abs(t.x), t.y, vec3(0, 0, 1));

	vec3 pos = vec3(t.x, 0, 0);
	col += plot(fbm_perlin(pos, 2., .5, .5), t.y, vec3(1, 0, 0));
	col += plot(fbm_simplex(pos, 2., .5, .5), t.y, vec3(0, 1, 0));
	//col += plot(fbm_worley(pos, 2., .5, .5), t.y, vec3(0, 0, 1));

#endif
// output
	fragColor = vec4 (col, 1);
}