#ifdef __cplusplus
#define _in(T) const T &
#define _inout(T) T &
#define _out(T) T &
#define _begin(type) type {
#define _end }
#define mul(a, b) (a) * (b)
#endif

#ifdef GL_ES
#define _in(T) const in T
#define _inout(T) inout T
#define _out(T) out T
#define _begin(type) type (
#define _end )
#define mul(a, b) (a) * (b)
precision mediump float;
#endif

#ifdef __HLSL
#define _in(T) const in T
#define _inout(T) inout T
#define _out(T) out T
#define _begin(type) {
#define _end }
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
#endif

#define PI 3.14159265359

// Shadertoy specific uniforms
#define u_res iResolution
#define u_time iGlobalTime
#define u_mouse iMouse

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
	float material_param;
	vec3 normal;
	vec3 origin;
};
#define max_dist 1e8
hit_t no_hit = _begin(hit_t)
	(max_dist + 1e1), // 'infinite' distance
	-1, // material id
	0., // material param
	vec3(0, 0, 0), // normal
	vec3(0, 0, 0) // origin
_end;

