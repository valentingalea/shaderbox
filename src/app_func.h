#include "def.h"

void mainImage(
	_out(vec4) fragColor,
#ifdef SHADERTOY
	vec2 fragCoord
#else
	_in(vec2) fragCoord
#endif
){
#define SCALE 4.
#define D (.025 * SCALE)

	vec2 t = fragCoord / u_res;
	t = (t* 2. - 1.) * SCALE;
	t.x += u_time * SCALE;
	
	float f = sin((t.x));
	
	f = smoothstep (f-D, f+D, t.y);
	f *= 1.- f;
	vec3 col = vec3 (f) * vec3 (5);
	
	fragColor = vec4 (col, 1);
}