// GLSL/C++ compatibility layer
#ifdef __cplusplus
#define _in(T) const T &
#define _inout(T) T &
#define _out(T) T &
#define _begin {
#define _end }
#else
#define _in(T) const in T
#define _inout(T) inout T
#define _out(T) out T
#define _begin (
#define _end )
#endif

#define PI 3.14159265359

// Shadertoy specific uniforms
#define u_res iResolution
#define u_time iGlobalTime
#define u_tex0 iChannel0

struct ray_t {
	vec3 origin;
	vec3 direction;
};
#define BIAS 1e-4 // small offset to avoid self-intersections

struct hit_t {
	float t;
	int material_id;
	float material_param;
	vec3 normal;
	vec3 origin;
};
#define max_dist 1e8 // TODO: precision issues
hit_t no_hit = hit_t _begin
	(max_dist + 1e1), -1, 1., vec3(0), vec3(0)
_end;

// camera system
vec3 eye, look_at;