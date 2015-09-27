// ----------------------------------------------------------------------------
// Noise functions
// ----------------------------------------------------------------------------
#include "../lib/ashima-noise/src/noise2d.glsl"

float fbm (
	_in(vec2) uv,
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
		total += amplitude * snoise (uv * frequency);
		amplitude *= gain;
		frequency *= lacunarity;
	}
	
	return total;
}