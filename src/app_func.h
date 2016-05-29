#include "def.h"

#include "../lib/ashima-noise/src/common.glsl"
#include "../lib/ashima-noise/src/classicnoise3d.glsl"
#include "../lib/ashima-noise/src/noise3d.glsl"
#include "../lib/ashima-noise/src/cellular3d.glsl"
#define noise(x) (snoise(x))

#include "fbm.h"
DECL_FBM_FUNC(fbm_perlin, 4, cnoise(p))
DECL_FBM_FUNC(fbm_simplex, 4, snoise(p))
DECL_FBM_FUNC(fbm_worley, 4, cellular(p).r)

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
	col += plot(fbm_worley(pos, 2., .5, .5), t.y, vec3(0, 0, 1));

// output
	fragColor = vec4 (col, 1);
}