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
	vec3 p = pos;
	
	t  = 0.51749673 * NOISE_FUNC(p); p *= 2.001;
	t += 0.25584929 * NOISE_FUNC(p); p *= 2.002;
	t += 0.12527603 * NOISE_FUNC(p); p *= 2.003;
	t += 0.06255931 * NOISE_FUNC(p);
	
	return t;
}

//TODO:
// - add rotation matrix for extra variation
// - domain distortion?