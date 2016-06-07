// ----------------------------------------------------------------------------
// Volumetric utilities
// ----------------------------------------------------------------------------

float isotropic_phase_func(float mu)
{
	return
	           1.
	/ //-------------------
	        4. * PI;
}

float rayleigh_phase_func(float mu)
{
	return
	        3. * (1. + mu*mu)
	/ //------------------------
	           (16. * PI);
}

// Henyey-Greenstein phase function factor [-1, 1]
// represents the average cosine of the scattered directions
// 0 is isotropic scattering
// > 1 is forward scattering, < 1 is backwards
//#define hg_g <external>

float henyey_greenstein_phase_func(float mu)
{
	return
	                     (1. - hg_g*hg_g)
	/ //---------------------------------------------
	     ((4. + PI) * pow(1. + hg_g*hg_g - 2.*hg_g*mu, 1.5));
}

float schlick_phase_func(float mu)
{
	// Schlick Phase Function factor
	// Pharr and  Humphreys [2004] equivalence to g from Henyey-Greenstein
#define shk_g (1.55*hg_g - 0.55 * (hg_g*hg_g*hg_g))

	return
	                  (1. - shk_g*shk_g)
	/ //-------------------------------------------
	     (4. * PI * (1. + shk_g*mu) * (1. + shk_g*mu));
}

struct volume_sampler_t {
	vec3 origin;			// start of ray
	vec3 pos;				// current pos of acccumulation ray
	float height;			// [0..1] within the volume
	float transmittance;	// (internal) energy loss by absorption & out-scattering
	vec3 radiance;			// mainly used as output color
	float alpha;
};

volume_sampler_t construct_volume(
	_in(vec3) origin
){
	volume_sampler_t v = _begin(volume_sampler_t)
		origin,
		origin,
		0.,
		1.,
		vec3(0, 0, 0),
		0.
	_end;
	return v;
}

