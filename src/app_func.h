#include "def.h"

float remap(
	_in(float) original_value,
	_in(float) original_min,
	_in(float) original_max,
	_in(float) new_min,
	_in(float) new_max
) {
	return new_min +
		(((original_value - original_min) /
		(original_max - original_min)) *
			(new_max - new_min));
}

#include "../lib/ashima-noise/src/common.glsl"
#include "../lib/ashima-noise/src/classicnoise3d.glsl"
#include "../lib/ashima-noise/src/noise3d.glsl"
#include "../lib/ashima-noise/src/cellular3d.glsl"
#define noise(x) (snoise(x))

#include "fbm.h"
DECL_FBM_FUNC(fbm_low_freq_perlin, 4, cnoise(p))
DECL_FBM_FUNC(fbm_low_freq_worley, 3, cellular(p).r)

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
	vec2 t = fragCoord.xy / u_res.xy;
#ifdef HLSL
	t.y = 1. - t.y;
#endif

	vec3 col = vec3(0, 0, 0);
#if 1 // 2D

	vec3 pos = vec3(t, 1);
	
	float p = fbm_low_freq_perlin(pos * 4., 2., .5, .5);
	float w = 1. - fbm_low_freq_worley(pos * 4., 4., .5, .5);
	float n = remap(p, -w, 1, 0, 1);

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
	float p = fbm_low_freq_perlin(pos * 4., 2., .5, .5);
	float w = 1. - fbm_low_freq_worley(pos * 4., 4., .5, .5);
	float n = remap(p, -w, 1, 0, 1);
	
	col += plot(n, t.y, vec3(0, 0, 1));
#endif

// output
	fragColor = vec4 (col, 1);
}