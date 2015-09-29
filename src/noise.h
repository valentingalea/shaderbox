// ----------------------------------------------------------------------------
// Noise functions
// ----------------------------------------------------------------------------
#include "../lib/ashima-noise/src/noise3d.glsl"

#define NOISE_FUNC snoise

float fbm_generic (
	_in(vec3) pos,
	_in(float) start_amplitude,
	_in(float) gain,
	_in(float) start_frequency,
	_in(float) lacunarity,
	_in(int) octaves
){
	float total = 0.;
	float amplitude = start_amplitude;
	float frequency = start_frequency;
	
	for (int i = 0; i < octaves; i++) {
		total += amplitude * NOISE_FUNC (pos * frequency);
		amplitude *= gain;
		frequency *= lacunarity;
	}
	
	return total;
}

float fbm_4(
	_in(vec3) pos	
){
	float t = 0.;
	float p = pos;
	
	t  = 0.5000 * NOISE_FUNC(p); p *= 2.001;
	t += 0.2500 * NOISE_FUNC(p); p *= 2.002;
	t += 0.1250 * NOISE_FUNC(p); p *= 2.003;
	t += 0.0625 * NOISE_FUNC(p);
	
	return t;
}

//TODO:
// - add rotation matrix for extra variation
// - for the unrolled on input more random amplitudes
// - domain distortion?
