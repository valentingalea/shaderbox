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
hit_t no_hit = hit_t _begin
	(max_dist + 1e1), // 'infinite' distance
	-1, // material id
	0., // material param
	vec3(0, 0, 0), // normal
	vec3(0, 0, 0) // origin
_end;

// camera system
vec3 eye, look_at;