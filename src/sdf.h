//
// Signed Distance Fields functions
//

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
