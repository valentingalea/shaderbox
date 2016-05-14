// ----------------------------------------------------------------------------
// Optics related utility functions
// ----------------------------------------------------------------------------

float fresnel_factor( // using Schlick’s approximation
	_in(float) n1,
	_in(float) n2,
	_in(float) VdotH // angle between View and half-dir between Light and View
){  
	float Rn = (n1 - n2) / (n1 + n2);
	float R0 = Rn * Rn; // reflection coeff. for light incoming parallel to the normal
	float F = 1. - VdotH;
	return R0 + (1. - R0) * (F * F * F * F * F);
}

#ifdef __cplusplus
vec3 reflect(
	_in(vec3) incident,
	_in(vec3) normal
){
	return incident - 2. * dot(normal, incident) * normal;
}

vec3 refract(
	_in(vec3) incident,
	_in(vec3) normal,
	_in(float) n
){
	float cosi = -dot(normal, incident);
	float sint2 = n * n * (1. - cosi * cosi);
	if (sint2 > 1.) {
		return reflect(incident, normal); // Total Internal Reflection - TODO: is this ok?
	}
	return n * incident + (n * cosi - sqrt(1. - sint2)) * normal;
}
#endif

