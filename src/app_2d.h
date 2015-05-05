#include "def.h"

sampler2D u_tex ("", sampler2D::Repeat);
#define u_res iResolution
#define u_time iGlobalTime

vec4 sample (_in(vec2) uv)
{
#if 1
	return texture (u_tex, uv);
#else
// TODO: some procedural shit
#endif
}

vec2 perturb (_in(vec2) uv)
{
	vec2 p = 2.*uv - 1.;
	
	float r = sqrt (dot (p, p));
	float a = atan (p.y, p.x);
	
	float s = cos (a) / r;
	float t = sin (a) / r;
	
	return vec2 (mod (s, 1.), mod (t, 1.));
}

void mainImage(_out(vec4) fragColor, _in(vec2) fragCoord)
{
	vec2 uv = fragCoord / u_res.xy;
	
	vec4 color = sample (
	perturb (uv) * 1.
	//+ u_time
	);
	
	fragColor = color;
}