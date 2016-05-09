// ----------------------------------------------------------------------------
// Fractal Brownian Motion
// depends on noise() primitive
// ----------------------------------------------------------------------------

#define DECL_FBM_FUNC(_name, _octaves, _gain) float _name(_in(vec3) pos, _in(float) lacunarity) { vec3 p = pos; float gain = _gain; float t = 0.; for (int i = 0; i < _octaves; i++) { t += gain * noise(p); p *= lacunarity; gain *= .5; } return t; }

#define DECL_TURB_FUNC(_name, _octaves, _gain) float _name(_in(vec3) pos, _in(float) lacunarity) { vec3 p = pos; float gain = _gain; float t = 0.; for (int i = 0; i < _octaves; i++) { t += gain * abs (noise(p)*2.-1.); p *= lacunarity; gain *= .5; } return t; }

DECL_FBM_FUNC(fbm, 4, .5)

float reference_fbm(
	_in(vec3) pos,
	_in(int) octaves,
	_in(float) lacunarity,
	_in(float) gain
){
	vec3 p = pos;
	float t = 0.;

	for (int i = 0; i < octaves; i++) {
		t += noise(p) * pow(lacunarity, -gain * i);
		p *= lacunarity;
	}

	return t;
}

float reference_multifractal(
	_in(vec3) pos,
	_in(int) octaves,
	_in(float) lacunarity,
	_in(float) gain,
	_in(float) offset
) {
	vec3 p = pos;
	float t = 1.;

	for (int i = 0; i < octaves; i++) {
		t *= (noise(p) + offset) * pow(lacunarity, -gain * i);
		p *= lacunarity;
	}

	return t;
}