// Modulo 289 without a division (only multiplications)
vec3 mod289(vec3 x)
{
	return x - floor(x * (1.0 / 289.0)) * 289.0;
}

vec4 mod289(vec4 x)
{
	return x - floor(x * (1.0 / 289.0)) * 289.0;
}

// Modulo 7 without a division
vec3 mod7(vec3 x)
{
	return x - floor(x * (1.0 / 7.0)) * 7.0;
}

vec4 mod7(vec4 x)
{
	return x - floor(x * (1.0 / 7.0)) * 7.0;
}

// Permutation polynomial: (34x^2 + x) mod 289
vec3 permute(vec3 x)
{
	return mod289((34.0 * x + 1.0) * x);
}

vec4 permute(vec4 x)
{
	return mod289(((x*34.0) + 1.0) * x);
}

vec4 taylorInvSqrt(vec4 r)
{
	return 1.79284291400159 - 0.85373472095314 * r;
}
