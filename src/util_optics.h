//
// Optics related utility functions
//

// Surfaces are more reflective at glancing angles
//
// how much reflection and how much transmission
// is given by Fresnel equations which can be
// Schlick's approximation
//    R(theta) = R0 + (1 - R0) (1 - cos theta)^5
// theta is incident angle
// cos(theta) = dot(N, V)
// R0 is the reflectance at normal incidence
// Can vary exponent for visual effect
//
// the direction of the refracted (aka transmitted)
// ray is given by the Snell eq given the
// indices of refraction of the mediums
//
// http://www.scratchapixel.com/old/lessons/3d-basic-lessons/lesson-14-interaction-light-matter/optics-reflection-and-refraction/
// http://graphics.stanford.edu/courses/cs148-10-summer/docs/2006--degreve--reflection_refraction.pdf
//    
// indices of refraction    
// gold ---- air ---- ice --- water --- beer --- alumim --- glass --- salt --- PET --- asphalt --- lead --- diamond ---- iron ---
// 0.470    1.000    1.309    1.333     1.345     1.390     1.500    1.516    1.575     1.645      2.010     2.420      2.950
//
float fresnel_factor(_in(float) n1, _in(float) n2, _in(float) VdotH)
{
	// using Schlick’s approximation    
	float Rn = (n1 - n2) / (n1 + n2);
	float R0 = Rn * Rn; // reflection coefficient for light incoming parallel to the normal
	float F = 1. - VdotH;
	return R0 + (1. - R0) * (F * F * F * F * F);
}

#ifdef __cplusplus
vec3 reflect(_in(vec3) incident, _in(vec3) normal)
{
	return incident - 2. * dot(normal, incident) * normal;
}

vec3 refract(_in(vec3) incident, _in(vec3) normal, _in(float) n)
{
	float cosi = -dot(normal, incident);
	float sint2 = n * n * (1. - cosi * cosi);
	if (sint2 > 1.) {
		return reflect(incident, normal); // Total Internal Reflection - TODO: is this ok?
	}
	return n * incident + (n * cosi - sqrt(1. - sint2)) * normal;
}
#endif
