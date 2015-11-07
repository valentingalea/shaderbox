#ifdef __cplusplus
#define _in(T) const T &
#define _inout(T) T &
#define _out(T) T &
#define _begin(type) type {
#define _end }
#define _mutable(T) T
#define _constant(T) const T
#define mul(a, b) (a) * (b)
#endif

#if defined(GL_ES) || defined(GL_SHADING_LANGUAGE_VERSION)
#define _in(T) const in T
#define _inout(T) inout T
#define _out(T) out T
#define _begin(type) type (
#define _end )
#define _mutable(T) T
#define _constant(T) const T
#define mul(a, b) (a) * (b)
precision mediump float;
#endif

#ifdef __HLSL
#define _in(T) const in T
#define _inout(T) inout T
#define _out(T) out T
#define _begin(type) {
#define _end }
#define _mutable(T) static T
#define _constant(T) static const T
#define vec2 float2
#define vec3 float3
#define vec4 float4
#define mat2 float2x2
#define mat3 float3x3
#define mat4 float4x4
#define mix lerp
#define fract frac
#define atan(y, x) atan2(x, y)
#define mod fmod
#pragma pack_matrix(row_major)
cbuffer uniforms : register(b0) {
	float2 u_res;
	float u_time;
	float2 u_mouse;
};
void mainImage(_out(float4) fragColor, _in(float2) fragCoord);
float4 main(float4 uv : SV_Position) : SV_Target{ float4 col; mainImage(col, uv.xy); return col; }
#endif

#if defined(__cplusplus) || defined(__SHADERTOY)
#define u_res iResolution
#define u_time iGlobalTime
#define u_mouse iMouse
#endif

#ifdef __GLSLSANDBOX
uniform float time;
uniform vec2 mouse;
uniform vec2 resolution;
#define u_res resolution
#define u_time time
#define u_mouse mouse
void mainImage(_out(vec4) fragColor, _in(vec2) fragCoord);
void main() { mainImage(gl_FragColor, gl_FragCoord.xy); }
#endif

#define PI 3.14159265359

struct ray_t {
	vec3 origin;
	vec3 direction;
};
#define BIAS 1e-4 // small offset to avoid self-intersections

struct sphere_t {
	vec3 origin;
	float radius;
	int material;
};

struct plane_t {
	vec3 direction;
	float distance;
	int material;
};

struct hit_t {
	float t;
	int material_id;
	vec3 normal;
	vec3 origin;
};
#define max_dist 1e8
hit_t no_hit = _begin(hit_t)
	float(max_dist + 1e1), // 'infinite' distance
	-1, // material id
	vec3(0., 0., 0.), // normal
	vec3(0., 0., 0.) // origin
_end;

