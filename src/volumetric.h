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
	vec3 origin; // start of ray
	vec3 pos; // current pos of acccumulation ray
	float height;

	float coeff_absorb;
	float T; // transmitance

	vec3 C; // color
	float alpha;
};

volume_sampler_t begin_volume(
	_in(vec3) origin,
	_in(float) coeff_absorb
){
	volume_sampler_t v = _begin(volume_sampler_t)
		origin, origin, 0.,
		coeff_absorb, 1.,
		vec3(0., 0., 0.), 0.
	_end;
	return v;
}

float illuminate_volume(
	_inout(volume_sampler_t) vol,
	_in(vec3) V,
	_in(vec3) L
);

void integrate_volume(
	_inout(volume_sampler_t) vol,
	_in(vec3) V,
	_in(vec3) L,
	_in(float) density,
	_in(float) dt
){
	// change in transmittance (follows Beer-Lambert law)
	float T_i = exp(-vol.coeff_absorb * density * dt);
	// Update accumulated transmittance
	vol.T *= T_i;
	// integrate output radiance (here essentially color)
	vol.C += vol.T * illuminate_volume(vol, V, L) * density * dt;
	// accumulate opacity
	vol.alpha += (1. - T_i) * (1. - vol.alpha);
}

