// ----------------------------------------------------------------------------
// Noise function by iq from https://www.shadertoy.com/view/ldl3Dl
// ----------------------------------------------------------------------------

vec3 hash_w(
	_in(vec3) x
){
#if 1
	vec3 xx = vec3(dot(x, vec3(127.1, 311.7, 74.7)),
		dot(x, vec3(269.5, 183.3, 246.1)),
		dot(x, vec3(113.5, 271.9, 124.6)));

	return fract(sin(xx)*43758.5453123);
#else
	return texture2D(iChannel0, (x.xy + vec2(3.0, 1.0)*x.z + 0.5) / 256.0, -100.0).xyz;
#endif
}

// returns closest, second closest, and cell id
vec3 noise_w(
	_in(vec3) x
){
	vec3 p = floor(x);
	vec3 f = fract(x);

	float id = 0.0;
	vec2 res = vec2(100.0, 100.0);
	for (int k = -1; k <= 1; k++)
		for (int j = -1; j <= 1; j++)
			for (int i = -1; i <= 1; i++)
			{
				vec3 b = vec3(float(i), float(j), float(k));
				vec3 r = vec3(b) - f + hash_w(p + b);
				float d = dot(r, r);

				if (d < res.x)
				{
					id = dot(p + b, vec3(1.0, 57.0, 113.0));
					res = vec2(d, res.x);
				}
				else if (d < res.y)
				{
					res.y = d;
				}
			}

	return vec3(sqrt(res), abs(id));
}
