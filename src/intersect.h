//
// Analytical surface-ray intersection routines
//

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

void intersect_sphere(_in(ray_t) ray, _in(sphere_t) sphere, _inout(hit_t) hit)
{
#if 1
	// geometrical solution
	// info: http://www.scratchapixel.com/old/lessons/3d-basic-lessons/lesson-7-intersecting-simple-shapes/ray-sphere-intersection/
	vec3 rc = sphere.origin - ray.origin;
	float radius2 = sphere.radius * sphere.radius;
	float tca = dot(rc, ray.direction);
	if (tca < 0.) return;
	float d2 = dot(rc, rc) - tca * tca;
	if (d2 > radius2) return;
	float thc = sqrt(radius2 - d2);
	float t0 = tca - thc;
	float t1 = tca + thc;

	if (t0 < 0.) t0 = t1;
	if (t0 < hit.t) {
#else // TODO: wrong for some reason... t gets weird values at intersection
	// analytical solution
	// based on combining the
	// sphere eq: (P - C)^2 = R^2
	// ray eq: P = O + t*D
	// into a quadratic eq: ax^2 + bx + c = 0
	// which can be solved by "completing the square" http://www.mathsisfun.com/algebra/completing-square.html
	// NOTE: be careful about "catastrophic cancellation" http://en.wikipedia.org/wiki/Loss_of_significance
	vec3 rc = ray.origin - sphere.origin;
	//	float a = D dot D -- which is 1 because D is normalised
	float b = 2.0 * dot(rc, ray.direction); // 2 * (O - C) dot D
	float c = dot(rc, rc) - sphere.radius * sphere.radius; // (O - C)^2 - R^2
	float discr = b * b - 4.0 * c;
	float t = -b - sqrt(abs(discr)); // use abs to avoid the check for < 0

	if (discr > 0.0 && t > 0.0 && t < hit.t) {
#endif
		vec3 impact = ray.origin + ray.direction * t0;

		hit.t = t0;
		hit.material_id = sphere.material;
		hit.material_param = 1.;
		hit.origin = impact;
		hit.normal = (impact - sphere.origin) / sphere.radius;
	}
}

// Plane is define by normal N and distance to origin P0 (which is on the plane itself)
// a plane eq is: (P - P0) dot N = 0
// which means that any line on the plane is perpendicular to the plane normal
// a ray eq: P = O + t*D
// substitution and solving for t gives:
// t = ((P0 - O) dot N) / (N dot D)
void intersect_plane(_in(ray_t) ray, _in(plane_t) p, _inout(hit_t) hit)
{
	float denom = dot(p.direction, ray.direction);
	if (denom > 1e-6)
	{
		float t = dot(vec3(p.distance) - ray.origin, p.direction) / denom;
		if (t >= 0.0 && t < hit.t)
		{
			vec3 impact = ray.origin + ray.direction * t;

			// checkboard pattern			
			//vec2 pattern = floor (impact.xz * 0.5);
			//float cb = mod (pattern.x + pattern.y, 2.0);

			hit.t = t;
			hit.material_id = p.material;
			hit.material_param = 1.; // cb; // Disabled for now
			hit.origin = impact;
			hit.normal = faceforward(p.direction, ray.direction, p.direction);
		}
	}
}
