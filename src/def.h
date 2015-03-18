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

struct ray_t {
	vec3 origin;
	vec3 direction;
};

struct hit_t {
	float t;
	int material_id;
	float material_param; // used by dynamic materials to store temp data
	vec3 normal;
	vec3 origin;
};

#define max_dist 999.0
hit_t no_hit = hit_t _begin
	(max_dist + 1.0), -1, 1., vec3(0), vec3(0)
_end;

#define BIAS 1e-4 // small offset to add to ray when retracing to avoid self-intersection
#define PI 3.14159265359

vec3 eye, look_at;