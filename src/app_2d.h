#include "def.h"
#include "noise.h"

sampler2D iChannel0 ("", sampler2D::Repeat);

#define u_tex iChannel0
#define u_res iResolution
#define u_time iGlobalTime

vec4 sample (_in(vec2) uv)
{
#if 1
	return texture (u_tex, uv);
#else
	float n = fbm (uv, 0.5, 0.5, 2., 2., 6);
	return vec4 (n, n, n, 1.);
#endif
}

vec2 perturb_road (_in(vec2) uv)
{
	vec2 p = 2.*uv - 1.;
	
	float s = p.x / abs(p.y);
	float t = 1. / abs(p.y);

	return  vec2 (s, t);
}

vec2 perturb_tunnel (_in(vec2) uv, _inout(float) r)
{
	vec2 p = 2.*uv - 1.;
	
	r = sqrt (dot (p, p));
	float a = atan (p.y, p.x) + u_time;
	
	float s = 1. / r + u_time;
	float t = 4. * (a / PI);

	return  vec2 (s, t);
}

float tent_filter (float t) // -1 to 1, peak at 0
{
	return max (1. - abs (t) , 0);
}

void mainImage(_out(vec4) fragColor, _in(vec2) fragCoord)
{
	vec2 uv = fragCoord / u_res.xy;
	
#if 1
	float r = 1./8.;
	uv.x += r * cos (u_time * 1.);
	uv.y += r * sin (u_time * 1.);
	
	float d;
	vec2 st = perturb_tunnel (uv, d);
	vec4 color = sample (st) * d;
#else
	vec2 st = perturb_road (uv);
	vec4 color = sample (st + vec2(0, u_time));

	color *= 1. - tent_filter (2.*uv.y - 1.);
#endif
		
	fragColor = color;
}