#include "def.h"

#if 0
#include "noise_iq.h"
#include "fbm.h"
#endif

#ifdef __cplusplus
sampler2D u_tex0 ("", sampler2D::Repeat);
#else
#define u_tex0 iChannel0
#endif

vec4 sample (
	_in(vec2) uv
){
#if 1
	return texture (u_tex0, uv);
#else
	float n = fbm (vec3(uv, 1.));
	return vec4 (n, n, n, 1.);
#endif
}

vec2 perturb_road (
	_in(vec2) uv,
	_in(float) time
){
	vec2 p = 2.*uv - 1.;
	
	float s = p.x / abs(p.y);
	float t = 1. / abs(p.y);

	return  vec2 (s, t - time);
}

vec2 perturb_tunnel (
	_in(vec2) uv,
	_in(float) time,
	_inout(float) r
){
	vec2 p = 2.*uv - 1.;
	
	r = sqrt (dot (p, p));
	float a = atan (p.y, p.x) + time;
	
	float s = 1. / r + time;
	float t = 4. * (a / PI);

	return  vec2 (s, t);
}

float tent_filter ( // -1 to 1, peak at 0
	_in(float) t
){
	return max (1. - abs (t) , 0);
}

void mainImage(
	_out(vec4) fragColor,
	_in(vec2) fragCoord
){
	vec2 uv = fragCoord / u_res.xy;
	
	vec2 st;
	vec4 color;
	float d = 1.;
	float t = mod (u_time, 16.);
	
	if (t < 4.) {
		st = perturb_tunnel (uv, u_time, d);
		color = sample (st);
		color *= d;
	}
	
	if (t > 4. && t < 8.) {
		st = perturb_tunnel (uv, 1., d);
		vec2 st2 = perturb_road (uv, 1.);
		color = sample (mix (st, st2, (t - 4.) / 4.));
		color *= d;	
	}
	
	if (t > 8. && t < 12.) {
		st = perturb_road (uv, u_time);
		color = sample (st);
	}

	if (t > 12.) {
		st = perturb_tunnel (uv, 1., d);
		vec2 st2 = perturb_road (uv, 1.);
		color = sample (mix (st2, st, (t - 12.) / 4.));
		color *= d;	
	}
	
	color *= 1. - tent_filter (2.*uv.y - 1.);

	// debug
	//color = vec4 (vec2 (fmod (st.x, 1.), fmod (st.y, 1.)),0, 1.);
		
	fragColor = color;
}