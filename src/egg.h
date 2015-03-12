#ifdef __cplusplus
#define SCREEN_WIDTH 100
#define SCREEN_HEIGHT 100
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
	float material_param;
	vec3 normal;
	vec3 origin;
};

#define max_dist 999.0
hit_t no_hit = hit_t _begin
	(max_dist + 1.0), -1, 1., vec3(0), vec3(0)
_end;

mat3 rotate_around_z(_in(float) angle_degrees);
mat3 rotate_around_y(_in(float) angle_degrees);
mat3 rotate_around_x(_in(float) angle_degrees);
vec3 corect_gamma(_in(vec3) color);
ray_t get_primary_ray(_in(vec3) cam_local_point, _in(vec3) cam_origin, _in(vec3) cam_look_at);
vec3 ik_solver(vec3 start, vec3 goal, float bone_length_1, float bone_length_2);

vec3 background(_in(ray_t) ray)
{
	return vec3(.1, .1, .7);
}

void setup_scene()
{
#define mat_debug 0
#define mat_egg 1
#define mat_bike 2
#define mat_ground 3
}

vec3 illuminate(_in(hit_t) hit)
{
#if 0 // debug: output the raymarching steps
	return vec3(hit.material_param);
#endif

	if (hit.material_id == mat_ground) return vec3(13. / 255., 104. / 255., 0. / 255.);
	if (hit.material_id == mat_egg) return vec3(0.95);
	if (hit.material_id == mat_bike) return vec3(.2);
	return vec3(1);
}

vec2 op_add(vec2 d1, vec2 d2) // union
{
	// minimum distance (preserving material info)
	return d1.x < d2.x ? d1 : d2;
}

float op_sub(float d1, float d2) // difference
{
	// intersection between first and
	// complement of the second field
	// aka the second 'carved out' from the first
	return max(d1, -d2);
}

float op_intersect(float d1, float d2)
{
	// what's common for both fields
	return max(d1, d2);
}

float op_blend(float a, float b, float k)
{
	// esentially liniar blend but more fancy
	// from http://iquilezles.org/www/articles/smin/smin.htm
	// NOTE: not true distance but estimate
	float h = clamp(0.5 + 0.5*(b - a) / k, 0.0, 1.0);
	return mix(b, a, h) - k*h*(1.0 - h);
}

float sd_plane(vec3 p, vec3 n, float d)
{
	// distance from point to plane
	// http://mathworld.wolfram.com/Point-PlaneDistance.html
	return dot(n, p) + d;
}

float sd_sphere(vec3 p, float s)
{
	// distance to center of sphere offset by the radius
	return length(p) - s;
}

float sd_box(vec3 p, vec3 b)
{
	// intersection of 3 axis aligned 'slabs'
	return max(abs(p.x) - b.x, max(abs(p.y) - b.y, abs(p.z) - b.z));
}

float sd_torus(vec3 p, float R, float r) // around z axis
{
	// projected circle of radius R on xy plane
	// combined with circle of radius r around z axis
	return length(vec2(length(p.xy) - R, p.z)) - r;
}

float sd_y_cylinder(vec3 p, float r, float h)
{
	// distance to the Y axis, offset (aka inflated) by the cylinder radius
	// then intersected with 2 cutting planes
	return max(length(p.xz) - r, abs(p.y) - h / 2.);
}

float sd_cylinder(vec3 P, vec3 P0, vec3 P1, float R)
{
	// distance to segment -- http://geomalgorithms.com/a02-_lines.html
	// then cut it with 2 planes at the ends
	// then offset it with radius    
	vec3 dir = normalize(P1 - P0);
	float dist = length(cross(dir, P - P0));
	float plane_1 = sd_plane(P, dir, length(P1));
	float plane_2 = sd_plane(P, -dir, -length(P0));
	return op_sub(op_sub(dist, plane_1), plane_2) - R;
}

vec2 sdf(_in(vec3) P)
{
	vec3 p = rotate_around_y(iGlobalTime * -50.0) * P -
		vec3(0, 0.5, 1.75);

	int material = mat_egg;

	float egg_y = 0.65;
#if 0
	float egg_m = sd_sphere(p - vec3(0, egg_y, 0), 0.475);
	float egg_b = sd_sphere(p - vec3(0, egg_y - 0.45, 0), 0.25);
	float egg_t = sd_sphere(p - vec3(0, egg_y + 0.45, 0), 0.25);
	float egg_1 = op_blend(egg_m, egg_b, .5);
	float egg_2 = op_blend(egg_1, egg_t, .5);
	vec2 egg = vec2(egg_2, material);
#else
    float s = 1.55;
    mat3 scale = mat3(
        s, 0, 0,
        0, 1, 0,
        0, 0, 1);
    mat3 iscale = mat3(
        1./s, 0, 0,
        0, 1./s, 0,
        0, 0, 1.);
    vec2 egg = vec2(
        sd_sphere(iscale * (scale * (p - vec3(0, egg_y, 0))), 0.475),
        material);
#endif

	vec3 wheel_pos = vec3(0, 1.2, 0);
	float pedal_radius = 0.3;
	float pedal_speed = 300.;
	float pedal_off = 0.2;

	mat3 rot_z = rotate_around_z(-iGlobalTime * pedal_speed);
	vec3 left_foot_pos = wheel_pos + rot_z * vec3(0, pedal_radius, pedal_off);

	rot_z = rotate_around_z(-iGlobalTime * pedal_speed);
	vec3 right_foot_pos = wheel_pos + rot_z * vec3(0, -pedal_radius, -pedal_off);

	vec3 side = vec3(0, 0, pedal_off);
	float femur = 0.8;
	float tibia = 0.75;
	float thick = .05;

	vec3 pelvis = vec3(0, 0., 0) + side;
	vec3 knee_l = ik_solver(pelvis, left_foot_pos, femur, tibia);
	vec2 left_leg_a = vec2(
		sd_cylinder(p + pelvis, vec3(0), knee_l - side, thick),
		material);
	vec2 left_leg_b = vec2(
		sd_cylinder(p + knee_l, vec3(0), left_foot_pos - knee_l, thick),
		material);

	pelvis = vec3(0, 0., 0) - side;
	vec3 knee_r = ik_solver(pelvis, right_foot_pos, femur, tibia);
	vec2 right_leg_a = vec2(
		sd_cylinder(p + pelvis, vec3(0), knee_r + side, thick),
		material);
	vec2 right_leg_b = vec2(
		sd_cylinder(p + knee_r, vec3(0), right_foot_pos - knee_r, thick),
		material);

	vec2 legs = op_add(
		vec2(op_blend(left_leg_a.x, left_leg_b.x, .01), material),
		op_add(right_leg_a, right_leg_b));

	vec3 left_toe = normalize(vec3(left_foot_pos.y - knee_l.y, knee_l.x - left_foot_pos.x, 0));
	vec2 left_foot = vec2(
		sd_cylinder(p + left_foot_pos, vec3(0), left_toe / 8., thick),
		material);

	vec3 right_toe = normalize(vec3(right_foot_pos.y - knee_r.y, knee_r.x - right_foot_pos.x, 0));
	vec2 right_foot = vec2(
		sd_cylinder(p + right_foot_pos, vec3(0), right_toe / 8., thick),
		material);

	vec2 feet = op_add(left_foot, right_foot);

	vec2 bike = vec2(
		sd_torus(p + wheel_pos, 1., .025),
		mat_bike);

	vec2 ground = vec2(
		sd_plane(P, vec3(0, 1, 0), wheel_pos.y + 0.5),
		mat_ground);

	return op_add(
		ground,
		op_add(legs, op_add(egg, op_add(feet, bike))));
}

vec3 sdf_normal(_in(vec3) p)
{
	float dt = 0.05;
	vec3 x = vec3(dt, 0, 0);
	vec3 y = vec3(0, dt, 0);
	vec3 z = vec3(0, 0, dt);
	return normalize(vec3(
		sdf(p + x).r - sdf(p - x).r,
		sdf(p + y).r - sdf(p - y).r,
		sdf(p + z).r - sdf(p - z).r
	));
}

#define EPSILON 0.001

float shadowmarch(_in(ray_t) ray) // from http://iquilezles.org/www/articles/rmshadows/rmshadows.htm
{
	const int steps = 20;
	const float end = 10.;
	const float penumbra_factor = 15.;
	const float darkest = 0.1;

	float t = 0.;
	float umbra = 1.;
	for (int i = 0; i < steps; i++) {
		vec3 p = ray.origin + ray.direction * t;
		vec2 d = sdf(p);

		if (t > end) break;
		if (d.x < EPSILON) {
			return darkest;
		}

		t += d.x;
		umbra = min(umbra, penumbra_factor * d.x / t);
	}

	return umbra;
}

vec3 raymarch(_in(ray_t) ray)
{
	const int steps = 75;
	const float end = 15.;

	float t = 0.;
	for (int i = 0; i < steps; i++) {
		vec3 p = ray.origin + ray.direction * t;
		vec2 d = sdf(p);

		if (t > end) break;
		if (d.x < EPSILON) {
			hit_t h = hit_t _begin
				t, // ray length at impact
				int(d.y), // material id
				float(i) / float(steps), // material custom param
				p, // point of impact
				vec3(0) // sdf_normal(p);
			_end;

			float s = 1.;
#if 1 // soft shadows
			if (int(d.y) == mat_ground) {
				vec3 sh_dir = vec3(0, 1, 1);
				ray_t sh_ray = ray_t _begin
					p + sh_dir * 0.05, sh_dir
				_end;
				s = shadowmarch(sh_ray);
			}
#endif

			return illuminate(h) * s;
		}

		t += d.x;
	}

	return background(ray);
}

void mainImage(_out(vec4) fragColor, _in(vec2) fragCoord)
{
	// assuming screen width is larger than height 
	vec2 aspect_ratio = vec2(iResolution.x / iResolution.y, 1);
	// field of view
	float fov = tan(radians(30.0));

	// antialising
#if 0
#define MSAA_PASSES 4
	float offset = 0.25;
	float ofst_x = offset * aspect_ratio.x;
	float ofst_y = offset;
	vec2 msaa[MSAA_PASSES];
	msaa[0] = vec2(-ofst_x, -ofst_y);
	msaa[1] = vec2(-ofst_x, +ofst_y);
	msaa[2] = vec2(+ofst_x, -ofst_y);
	msaa[3] = vec2(+ofst_x, +ofst_y);
#else
#define MSAA_PASSES 1
	vec2 msaa[MSAA_PASSES];
	msaa[0] = vec2(0.5);
#endif

	vec3 eye = vec3(1, 2, 5);
	vec3 look_at = vec3(0);

	vec3 color = vec3(0);

	for (int i = 0; i < MSAA_PASSES; i++) {
		vec2 point_ndc = (fragCoord.xy + msaa[i]) / iResolution.xy;
		vec3 point_cam = vec3((2.0 * point_ndc - 1.0) * aspect_ratio * fov, -1.0);

		ray_t ray = get_primary_ray(point_cam, eye, look_at);

		color += raymarch(ray) / float(MSAA_PASSES);
	}

	fragColor = vec4(corect_gamma(color), 1);
}

mat3 rotate_around_z(_in(float) angle_degrees)
{
	float angle = radians(angle_degrees);
	float _sin = sin(angle);
	float _cos = cos(angle);
	return mat3(_cos, -_sin, 0, _sin, _cos, 0, 0, 0, 1);
}

mat3 rotate_around_y(_in(float) angle_degrees)
{
	float angle = radians(angle_degrees);
	float _sin = sin(angle);
	float _cos = cos(angle);
	return mat3(_cos, 0, _sin, 0, 1, 0, -_sin, 0, _cos);
}

mat3 rotate_around_x(_in(float) angle_degrees)
{
	float angle = radians(angle_degrees);
	float _sin = sin(angle);
	float _cos = cos(angle);
	return mat3(1, 0, 0, 0, _cos, -_sin, 0, _sin, _cos);
}

vec3 corect_gamma(_in(vec3) color)
{
	float gamma = 1.0 / 2.25;
	return vec3(pow(color.r, gamma), pow(color.g, gamma), pow(color.b, gamma));
}

ray_t get_primary_ray(_in(vec3) cam_local_point, _in(vec3) cam_origin, _in(vec3) cam_look_at)
{
	vec3 fwd = normalize(cam_look_at - cam_origin);
	vec3 up = vec3(0, 1, 0);
	vec3 right = cross(up, fwd);
	up = cross(fwd, right);

	return ray_t _begin
		cam_origin,
		normalize(fwd + up * cam_local_point.y + right * cam_local_point.x)
	_end;
}

vec3 ik_2_bone_centered_solver(vec3 goal, float L1, float L2)
{
#if 0 // from https://www.shadertoy.com/view/ldlGR7
	vec3 q = goal*(0.5 + 0.5*(L1*L1 - L2*L2) / dot(goal, goal));

	float s = L1*L1 - dot(q, q);
	s = max(s, 0.0);
	q += sqrt(s)*normalize(cross(goal, vec3(0, 0, 1)));

	return q;
#else // naive version with law of cosines
	float G = length(goal);

	// tetha is the angle between bone1 and goal direction
	// get it from law of cosines applied to the
	// triangle with sides: bone1, bone2, pivot_of_bone1<->goal
	float cos_theta = (L1*L1 + G*G - L2*L2) / (2.*L1*G);

	// sin^2 + cos^2 = 1 (Pythagoras in unit circle)
	float sin_theta = sqrt(1. - cos_theta * cos_theta);

	// rotation matrix by theta amount around the axis
	// perpendicular to the plane created by bone1 and bone2
	mat3 rot = mat3(
		cos_theta, -sin_theta, 0,
		sin_theta, cos_theta, 0,
		0, 0, 1.
		);

	// get the end of bone1 aka the pivot of bone2
	// by getting a vector from the goal direction
	// and rotating along with the newly found theta angle
	return rot * (normalize(goal) * L1);
#endif
}

vec3 ik_solver(vec3 start, vec3 goal, float bone_length_1, float bone_length_2)
{
	return start + ik_2_bone_centered_solver(
		goal - start, bone_length_1, bone_length_2);
}
