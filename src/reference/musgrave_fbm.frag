/////////////////////////////////////////////////////////////////////////////////////////////
//
// Musgrave's noise functions for terrain rendering 
// https://engineering.purdue.edu/~ebertd/texture/1stEdition/musgrave/musgrave.c
// Originaly published in "Texturing & Modeling: A Procedural Approach"
// by David S. Ebert, F. Kenton Musgrave, Darwyn Peachey, Ken Perlin, Steven Worley
//
// Another Musgrave Paper
// http://www.classes.cs.uchicago.edu/archive/2014/winter/23700-1/project_4_and_5/MusgraveTerrain00.pdf
//
// Use Ashima Simplex noise Implementation
// https://github.com/ashima/webgl-noise
//
// License: Distributed under the MIT License.
//
// Include Ashima port of Musgrave code:
// - fBmA(),
// - multifractalA(),
// - heteroTerrainA()
//
// Ashima port of multifractalA() and heteroTerrainA() differs from original code.
// It corrects some bugs, adds a new parameter for lacunarity and handles lacunarity
// differently than Musgrave (in a more standard way IMHO).
// I kept this modification in my personal port of the original code. Although I tried
// to be closer to the original.
//
// My port of Musgrave code:
// - fBm(),
// - multifractal(),
// - heteroTerrain(),
// - hybridMultiFractal(),
// - ridgedMultifractal()
//
// Displays from left to right:
// - fBm()
// - multifractal()
// - heteroTerrain()
// - hybridMultiFractal()
// - ridgedMultifractal()
//
// Top line: function result scaled to 0..1 for grayscale output
// Bottom line: function result map to Yellow for positive values and Cyan for negative values.
//				(black means noise == 0. at this point)
//

/////////////////////////////////////////////////////////////////////////////////////////////
// Ashima Simplex Noise
//
// Description : Array and textureless GLSL 2D simplex noise function.
//      Author : Ian McEwan, Ashima Arts.
//  Maintainer : ijm
//     Lastmod : 20110822 (ijm)
//     License : Copyright (C) 2011 Ashima Arts. All rights reserved.
//               Distributed under the MIT License. See LICENSE file.
//               https://github.com/ashima/webgl-noise
// 

vec3 mod289(vec3 x) {
  return x - floor(x * (1.0 / 289.0)) * 289.0;
}

vec2 mod289(vec2 x) {
  return x - floor(x * (1.0 / 289.0)) * 289.0;
}

vec3 permute(vec3 x) {
  return mod289(((x*34.0)+1.0)*x);
}

float snoise(vec2 v)
  {
  const vec4 C = vec4(0.211324865405187,  // (3.0-sqrt(3.0))/6.0
                      0.366025403784439,  // 0.5*(sqrt(3.0)-1.0)
                     -0.577350269189626,  // -1.0 + 2.0 * C.x
                      0.024390243902439); // 1.0 / 41.0
// First corner
  vec2 i  = floor(v + dot(v, C.yy) );
  vec2 x0 = v -   i + dot(i, C.xx);

// Other corners
  vec2 i1;
  //i1.x = step( x0.y, x0.x ); // x0.x > x0.y ? 1.0 : 0.0
  //i1.y = 1.0 - i1.x;
  i1 = (x0.x > x0.y) ? vec2(1.0, 0.0) : vec2(0.0, 1.0);
  // x0 = x0 - 0.0 + 0.0 * C.xx ;
  // x1 = x0 - i1 + 1.0 * C.xx ;
  // x2 = x0 - 1.0 + 2.0 * C.xx ;
  vec4 x12 = x0.xyxy + C.xxzz;
  x12.xy -= i1;

// Permutations
  i = mod289(i); // Avoid truncation effects in permutation
  vec3 p = permute( permute( i.y + vec3(0.0, i1.y, 1.0 ))
		+ i.x + vec3(0.0, i1.x, 1.0 ));

  vec3 m = max(0.5 - vec3(dot(x0,x0), dot(x12.xy,x12.xy), dot(x12.zw,x12.zw)), 0.0);
  m = m*m ;
  m = m*m ;

// Gradients: 41 points uniformly over a line, mapped onto a diamond.
// The ring size 17*17 = 289 is close to a multiple of 41 (41*7 = 287)

  vec3 x = 2.0 * fract(p * C.www) - 1.0;
  vec3 h = abs(x) - 0.5;
  vec3 ox = floor(x + 0.5);
  vec3 a0 = x - ox;

// Normalise gradients implicitly by scaling m
// Approximation of: m *= inversesqrt( a0*a0 + h*h );
  m *= 1.79284291400159 - 0.85373472095314 * ( a0*a0 + h*h );

// Compute final noise value at P
  vec3 g;
  g.x  = a0.x  * x0.x  + h.x  * x0.y;
  g.yz = a0.yz * x12.xz + h.yz * x12.yw;
  return 130.0 * dot(m, g);
}

///////////////////////////////////////////////////////////////////////////////////////////
// K.Musgrave Procedural Noises Collections
// Ashima Port
//

/*
 * Procedural fBm evaluated at “point”.
 * 
 * Parameters:
 * “H” is the fractal increment parameter
 * “lacunarity” is the gap between successive frequencies
 * “octaves” is the number of frequencies in the fBm
 *
 * Ebert, D., F. K. Musgrave, D. Peachey, K. Perlin, and S. Worley. 2003. Texturing and modeling: A procedural approach, 437. Third Edition. San Francisco: Morgan Kaufmann.
*/
float fBmA(vec2 point, float H, float lacunarity, float frequency, float octaves) {
  float value = 0.0;
  float remainder = 0.0;
  float pwrHL = pow(lacunarity, -H);
  float pwr = 1.0; 

  /* inner loop of fractal construction */
  for (int i=0; i<65535; i++) {
    value += snoise(point * frequency) * pwr;
    pwr *= pwrHL;
    point *= lacunarity;

    if (i==int(octaves)-1) break;
  }

  remainder = octaves - floor(octaves);
  if (remainder != 0.0) {
    value += remainder * snoise(point * frequency) * pwr;
  }

  return value;
}

/*
 * Procedural multifractal evaluated at “point.”
 * 
 * Parameters:
 * “H” determines the fractal dimension of the roughest areas 
 * “lacunarity” is the gap between successive frequencies 
 * “octaves” is the number of frequencies in the fBm
 * “offset” is the zero offset, which determines multifractality
 *
 * Ebert, D., F. K. Musgrave, D. Peachey, K. Perlin, and S. Worley. 2003. Texturing and modeling: A procedural approach, 440. Third Edition. San Francisco: Morgan Kaufmann.
*/
float multifractalA(vec2 point, float H, float lacunarity, float frequency, float octaves, float offset) {
  float value = 1.0;
  float rmd = 0.0;
  float pwHL = pow(lacunarity, -H);
  float pwr = 1.0;
  
  /* inner loop of fractal construction */
  for (int i=0; i<65535; i++) {
    value *= pwr * snoise(point*frequency) + offset;
    pwr *= pwHL;
    point *= lacunarity;

    if (i==int(octaves)-1) break;
  }

  rmd = octaves - floor(octaves);
  if (rmd != 0.0) value += (rmd * snoise(point*frequency) * pwr);

  return value;
}

/*
 * Heterogeneous procedural terrain function: stats by altitude method. 
 * Evaluated at “point”; returns value stored in “value”.
 * 
 * Parameters:
 * “H” determines the fractal dimension of the roughest areas 
 * “lacunarity” is the gap between successive frequencies 
 * “octaves” is the number of frequencies in the fBm
 * “offset” raises the terrain from “sea level”
 *
 * Ebert, D., F. K. Musgrave, D. Peachey, K. Perlin, and S. Worley. 2003. Texturing and modeling: A procedural approach, 500. Third Edition. San Francisco: Morgan Kaufmann.
*/
float heteroTerrainA(vec2 point, float H, float lacunarity, float frequency, float octaves, float offset) {
 float value, increment, remainder;
  float pwrHL = pow(lacunarity, -H);
  float pwr = pwrHL; /* starts with i=1 instead of 0 */

  value = offset + snoise(point * frequency);      
  point *= lacunarity;

  for (int i=1; i<65535; i++) {
    increment = (snoise(point * frequency) + offset) * pwr * value;
    // frequency *= lacunarity;
    value += increment;
    point *= lacunarity;

    if (i==int(octaves)) break;
  }

  /* take care of remainder in 'octaves'  */
  remainder = mod(octaves, floor(octaves));

  if (remainder != 0.0) {
    increment = (snoise(point * frequency) + offset) * pwr * value;
    value += remainder * increment;
  }

  return value;  
}

///////////////////////////////////////////////////////////////////////////////////////////
// K.Musgrave Procedural Noises Collections
// XBE Port
//

/*
 * Procedural fBm evaluated at “point”.
 * 
 * Parameters:
 * “H” is the fractal increment parameter
 * “lacunarity” is the gap between successive frequencies
 * “octaves” is the number of frequencies in the fBm
 *
 * Ebert, D., F. K. Musgrave, D. Peachey, K. Perlin, and S. Worley. 2003. Texturing and modeling: A procedural approach, 437. Third Edition. San Francisco: Morgan Kaufmann.
*/
float fBm(vec2 point, float H, float lacunarity, float frequency, float octaves)
{
	float value = 0.0;
	float rmd = 0.0;
	float pwHL = pow(lacunarity, -H);
	float pwr = pwHL; 

	for (int i=0; i<65535; i++)
	{
		value += snoise(point * frequency) * pwr;
		point *= lacunarity;
		pwr *= pwHL;
		if (i==int(octaves)-1) break;
	}

	rmd = octaves - floor(octaves);
	if (rmd != 0.0) value += rmd * snoise(point * frequency) * pwr;

	return value;
}

/*
 * Procedural multifractal evaluated at “point.”
 * 
 * Parameters:
 * “H” determines the fractal dimension of the roughest areas 
 * “lacunarity” is the gap between successive frequencies 
 * “octaves” is the number of frequencies in the fBm
 * “offset” is the zero offset, which determines multifractality
 *
 * Ebert, D., F. K. Musgrave, D. Peachey, K. Perlin, and S. Worley. 2003. Texturing and modeling: A procedural approach, 440. Third Edition. San Francisco: Morgan Kaufmann.
*/
float multifractal(vec2 point, float H, float lacunarity, float frequency, float octaves, float offset)
{
	float value = 1.0;
	float rmd = 0.0;
	float pwHL = pow(lacunarity, -H);
	float pwr = pwHL;

	for (int i=0; i<65535; i++)
	{
		value *= pwr*snoise(point*frequency) + offset;
		point *= lacunarity;
		pwr *= pwHL;
		if (i==int(octaves)-1) break;
	}

	rmd = octaves - floor(octaves);
	if (rmd != 0.0) value += (rmd * snoise(point*frequency) * pwr);

	return value;
}

/*
 * Heterogeneous procedural terrain function: stats by altitude method. 
 * Evaluated at “point”; returns value stored in “value”.
 * 
 * Parameters:
 * “H” determines the fractal dimension of the roughest areas 
 * “lacunarity” is the gap between successive frequencies 
 * “octaves” is the number of frequencies in the fBm
 * “offset” raises the terrain from “sea level”
 *
 * Ebert, D., F. K. Musgrave, D. Peachey, K. Perlin, and S. Worley. 2003. Texturing and modeling: A procedural approach, 500. Third Edition. San Francisco: Morgan Kaufmann.
*/
float heteroTerrain(vec2 point, float H, float lacunarity, float frequency, float octaves, float offset)
{
	float value = 1.;
	float increment = 0.;
	float rmd = 0.;
	float pwHL = pow(lacunarity, -H);
	float pwr = pwHL;

	value = pwr*(offset + snoise(point * frequency));
	point *= lacunarity;
	pwr *= pwHL;

	for (int i=1; i<65535; i++)
	{
		increment = (snoise(point * frequency) + offset) * pwr * value;
		value += increment;
		point *= lacunarity;
		pwr *= pwHL;
		if (i==int(octaves)) break;
	}

	rmd = mod(octaves, floor(octaves));
	if (rmd != 0.0) value += rmd * ((snoise(point * frequency) + offset) * pwr * value);

	return value;  
}

/* Hybrid additive/multiplicative multifractal terrain model.
 *
 * Copyright 1994 F. Kenton Musgrave 
 *
 * Some good parameter values to start with:
 *
 *      H:           0.25
 *      offset:      0.7
 */
float hybridMultiFractal(vec2 point, float H, float lacunarity, float frequency, float octaves, float offset)
{
	float value = 1.0;
	float signal = 0.0;
	float rmd = 0.0;
	float pwHL = pow(lacunarity, -H);
	float pwr = pwHL;
	float weight = 0.;

	/* get first octave of function */
	value = pwr*(snoise(point * frequency)+offset);
	weight = value;
	point *= lacunarity;
	pwr *= pwHL;

	/* spectral construction inner loop, where the fractal is built */
	for (int i=1; i<65535; i++)
	{
		weight = weight>1. ? 1. : weight;
		signal = pwr * (snoise(point*frequency) + offset);
		value += weight*signal;
		weight *= signal;
		pwr *= pwHL;
		point *= lacunarity;
		if (i==int(octaves)-1) break;
	}

	/* take care of remainder in ``octaves''  */
	rmd = octaves - floor(octaves);
	if (rmd != 0.0) value += (rmd * snoise(point*frequency) * pwr);

	return value;
}

/* Ridged multifractal terrain model.
 *
 * Copyright 1994 F. Kenton Musgrave 
 *
 * Some good parameter values to start with:
 *
 *      H:           1.0
 *      offset:      1.0
 *      gain:        2.0
 */
float ridgedMultiFractal(vec2 point, float H, float lacunarity, float frequency, float octaves, float offset, float gain)
{
	float value = 1.0;
	float signal = 0.0;
	float pwHL = pow(lacunarity, -H);
	float pwr = pwHL;
	float weight = 0.;

	/* get first octave of function */
	signal = snoise(point * frequency);
	signal = offset-abs(signal);
	signal *= signal;
	value = signal * pwr;
	weight = 1.0;
	pwr *= pwHL;

	/* spectral construction inner loop, where the fractal is built */
	for (int i=1; i<65535; i++)
	{
		point *= lacunarity;
		weight = clamp(signal*gain, 0.,1.);
		signal = snoise(point * frequency);
		signal = offset-abs(signal);
		signal *= signal;
		signal *= weight;
		value += signal * pwr;
		pwr *= pwHL;
		if (i==int(octaves)-1) break;
	}

	return value;
}

/////////////////////////////////////////////////////////////////////////////////////////////
//

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
	vec2 uv = fragCoord.xy / iResolution.xy;
	vec2 q = 2.*uv-1.;
	q.y = q.y<0.?1.+q.y:q.y;
	q *= 2.;
	// Noise functions
	float n = -9.;
	if (uv.x < 0.198)
		n = fBm(q, 0.5, 2., 2., 6.);
	if ((uv.x > 0.2) && ( uv.x < 0.398 ))
	{
		float o = 0.65;
		n = multifractal(q, 0.5, 2., 2., 6., o)/o;
	}
	if ((uv.x > 0.4) && ( uv.x < 0.598 ))
		n = heteroTerrain(q, 0.5, 2., 2., 6., 0.05);
	if ((uv.x > 0.6) && ( uv.x < 0.798 ))
		n = hybridMultiFractal(q, 0.25, 3., 2., 6., 0.05);
	if (uv.x>0.8)
		n = ridgedMultiFractal(q, 0.25, 2., 2., 6., .9, 1.5);
	// Display
	if (n<-8.9)
		discard;
	vec3 col = vec3(0.0);
	if (uv.y>0.502)
	{
		n = 0.5+0.5*n;
		col = vec3(clamp(n,0.,1.));
	}
	else if (uv.y < 0.5)
	{
		col = vec3(clamp(n,0.,1.), clamp(n,0.,1.)+clamp(-n,0.,1.), clamp(-n,0.,1.) );
		col = sqrt(col);
	}
	fragColor = vec4(col,1.0);
}