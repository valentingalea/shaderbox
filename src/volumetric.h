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
