// ----------------------------------------------------------------------------
// Volumetric utilities
// ----------------------------------------------------------------------------

// Notes taken from
// Production Volume Rendering by Magnus Wrenninge
//
// *participating media*
// for integration purposes we usually consider
// a differential cylinder filled with a participating medium
// TODO: drawing
// TODO: notations
// dL(p, ω) = Li(p,ω) − Lo(p,ω)
//          = emission + scattering_in − scattering_out − absorption
//
// *absorption*
// The unit of the absorption coefficient σa is a 
// reciprocal distance (i.e. m-1), which essentially
// describes the probability that a ray of light
// is absorbed as it travels through one length
// unit of the medium.
// Being a reciprocal means it can be infinitely large.
//
// Lo = Li + dL
// dL = -σa  Li  dt
//
// *emission*
// <not used here>
//
// *scaterring*
// The scattering property describes the
// likelihood that a ray of light traveling
// through the medium will bounce and reflect
// to another direction.
//
// this interaction accounts for both
// inscattering and out-scattering.
//
// Lo = Li  +  dLin  + dLout
// dLin  =  σs  phase(ω,  ω') Light(p,  ω')
// dLout  = -σs  Li  dt
//
// *extinction*
// When treated together, out-scattering and
// absorption are referred to as extinction or
// attenuation, and it uses the symbol σe.
//
// *phase functions*
// A black box indicating how light how much
// light traveling through a medium in
// direction  ω  will, upon scattering,
// reflect to direction  ω', i.e. 
// probability = p(ω,  ω').
//
// Similar to BRDF except
// * it usually considers a single angle
//   (the phase angle between 2 directions)
// * integrates to 1 over the entire sphere of directions
//
// *tranmission*
// The relationship L(s) / L(0) is called transmittance
// and it is the definition of how much light
// can pass between two points in the medium
//
// The incoming and outgoing radiance measures
// are directly related by transmittance
// and the relationship is called Beer’s law
// where T is the transmittance and 
// τ is the optical thickness measure
//
// T = exp (−τ)
// τ = (σs + σa)s

float rayleigh_phase_func(float mu)
{
	return
	        3. * (1. + mu*mu)
	/ //------------------------
	           (16. * PI);
}

float henyey_greenstein_phase_func(float mu)
{
	// Henyey-Greenstein phase function factor [-1, 1]
	// represents the average cosine of the scattered directions
	// 0 is isotropic scattering
	// > 1 is forward scattering, < 1 is backwards
	const float g = 0.76;

	return
	                     (1. - g*g)
	/ //---------------------------------------------
	     ((4. + PI) * pow(1. + g*g - 2.*g*mu, 1.5));
}

float schlick_phase_func(float mu)
{
	// Schlick Phase Function factor
	// Pharr and  Humphreys [2004] equivalence to g from Henyey-Greenstein
	const float g = 0.76;
	const float k = 1.55*g - 0.55 * (g*g*g);

	return
	                  (1. - k*k)
	/ //-------------------------------------------
	     (4. * PI * (1. + k*mu) * (1. + k*mu));
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

