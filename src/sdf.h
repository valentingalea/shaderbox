// ----------------------------------------------------------------------------
// Signed Distance Fields functions
// ----------------------------------------------------------------------------

vec2 op_add( // union
	_in(vec2) d1,
	_in(vec2) d2
){
	// minimum distance (preserving material info)
	return d1.x < d2.x ? d1 : d2;
}

float op_sub( // difference
	_in(float) d1,
	_in(float) d2
){
	// intersection between first and
	// complement of the second field
	// aka the second 'carved out' from the first
	return max(d1, -d2);
}

float op_intersect( // intersection
	_in(float) d1,
	_in(float) d2
){
	// what's common for both fields
	return max(d1, d2);
}

float op_blend(
	_in(float) a,
	_in(float) b,
	_in(float) k // factor of smoothing
){
	// from http://iquilezles.org/www/articles/smin/smin.htm
	// NOTE: not true distance but estimate
	float h = clamp(0.5 + 0.5*(b - a) / k, 0.0, 1.0);
	return mix(b, a, h) - k*h*(1.0 - h);
}

float sd_plane(
	_in(vec3) p,
	_in(vec3) n, // normal
	_in(float) d // distance
){
	// distance from point to plane
	// http://mathworld.wolfram.com/Point-PlaneDistance.html
	return dot(n, p) + d;
}

float sd_sphere(
	_in(vec3) p,
	_in(float) r
){
	// distance to center of sphere offset by the radius
	return length(p) - r;
}

float sd_box(
	_in(vec3) p,
	_in(vec3) b // dimensions of box
){
	// intersection of 3 axis aligned 'slabs'
	return max(abs(p.x) - b.x, max(abs(p.y) - b.y, abs(p.z) - b.z));
}

float sd_torus( // around Z axis
	_in(vec3) p,
	_in(float) R, // 'donut' radius
	_in(float) r  // thickness
){
	// projected circle of radius R on xy plane
	// combined with circle of radius r around z axis
	return length(vec2(length(p.xy) - R, p.z)) - r;
}

float sd_y_cylinder(
	_in(vec3) p,
	_in(float) r, // radius
	_in(float) h  // height
){
	// distance to the Y axis, offset (aka inflated) by the cylinder radius
	// then intersected with 2 cutting planes
	return max(length(p.xz) - r, abs(p.y) - h / 2.);
}

float sd_cylinder(
	_in(vec3) P,
	_in(vec3) P0, // start point
	_in(vec3) P1, // end point
	_in(float) R  // thickness
){
	// distance to segment -- http://geomalgorithms.com/a02-_lines.html
	// then cut it with 2 planes at the ends
	// then offset it with radius    
	vec3 dir = normalize(P1 - P0);
	float dist = length(cross(dir, P - P0));
	float plane_1 = sd_plane(P, dir, length(P1));
	float plane_2 = sd_plane(P, -dir, -length(P0));
	return op_sub(op_sub(dist, plane_1), plane_2) - R;
}

// 3D Bezier curved cylinder
// original by http://research.microsoft.com/en-us/um/people/hoppe/ravg.pdf
// adapted by iq https://www.shadertoy.com/view/ldj3Wh
float det(
	_in(vec2) a,
	_in(vec2) b
){
	return a.x*b.y - b.x*a.y;
}
vec3 sd_bezier_get_closest(
	_in(vec2) b0,
	_in(vec2) b1,
	_in(vec2) b2
){
	float a = det(b0, b2);
	float b = 2.0*det(b1, b0);
	float d = 2.0*det(b2, b1);
	float f = b*d - a*a;
	vec2  d21 = b2 - b1;
	vec2  d10 = b1 - b0;
	vec2  d20 = b2 - b0;
	vec2  gf = 2.0*(b*d21 + d*d10 + a*d20); gf = vec2(gf.y, -gf.x);
	vec2  pp = -f*gf / dot(gf, gf);
	vec2  d0p = b0 - pp;
	float ap = det(d0p, d20);
	float bp = 2.0*det(d10, d0p);
	float t = clamp((ap + bp) / (2.0*a + b + d), 0.0, 1.0);
	return vec3(mix(mix(b0, b1, t), mix(b1, b2, t), t), t);
}
vec2 sd_bezier(
	_in(vec3) a, // start
	_in(vec3) b, // knot (control point)
	_in(vec3) c, // end
	_in(vec3) p, 
	_in(float) thickness
){
	vec3 w = normalize(cross(c - b, a - b));
	vec3 u = normalize(c - b);
	vec3 v = normalize(cross(w, u));

	vec2 a2 = vec2(dot(a - b, u), dot(a - b, v));
	vec2 b2 = vec2(0.0);
	vec2 c2 = vec2(dot(c - b, u), dot(c - b, v));
	vec3 p3 = vec3(dot(p - b, u), dot(p - b, v), dot(p - b, w));

	vec3 cp = sd_bezier_get_closest(a2 - p3.xy, b2 - p3.xy, c2 - p3.xy);

	return vec2(0.85*(sqrt(dot(cp.xy, cp.xy) + p3.z*p3.z) - thickness), cp.z);
}