// ----------------------------------------------------------------------------
// Fractal Brownian Motion
// ----------------------------------------------------------------------------

float fbm(
	_in(vec3) pos,
	_in(float) lacunarity
){
	vec3 p = pos;
	float
	t  = 0.51749673 * noise(p); p *= lacunarity;
	t += 0.25584929 * noise(p); p *= lacunarity;
	t += 0.12527603 * noise(p); p *= lacunarity;
	t += 0.06255931 * noise(p);
	
	return t;
}