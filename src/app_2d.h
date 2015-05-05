#include "def.h"

sampler2D iChannel0 ("", sampler2D::Repeat);

#define u_tex iChannel0
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
	
//	float s = p.x*cos(2.*r) - p.y*sin(2.*r);
//	float t = p.x*cos(2.*r) - p.y*sin(2.*r);
//  float s = 1./(r+0.5+0.5*sin(5.*a));
//	float t = a*3./PI;
	float s = p.x / abs(p.y);
	float t = 1. / abs(p.y);

	return vec2 (mod (s, 1.), mod (t, 1.));
}

void mainImage(_out(vec4) fragColor, _in(vec2) fragCoord)
{
	vec2 uv = fragCoord / u_res.xy;
	
	vec4 color = sample (
		perturb (uv)
		+ vec2(0, u_time)
	);
	
	fragColor = color;
}