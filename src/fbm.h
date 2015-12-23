// ----------------------------------------------------------------------------
// Fractal Brownian Motion
// ----------------------------------------------------------------------------

// general algorithm:
//	for (int i = 0; i < octaves; i++) {
//		total += amplitude * noise (pos * frequency);
//		amplitude *= gain;
//		frequency *= lacunarity;
//	}
	
float fbm(
	_in(vec3) pos	
){
	vec3 p = pos;
	float
	t  = 0.51749673 * noise(p); p *= 3.003;
	t += 0.25584929 * noise(p); p *= 3.001;
	t += 0.12527603 * noise(p); p *= 3.002;
	t += 0.06255931 * noise(p);
	
	return t;
}