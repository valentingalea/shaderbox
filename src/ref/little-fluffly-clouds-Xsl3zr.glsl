// little fluffy clouds
// @simesgreen

const int _VolumeSteps = 64;
const float _StepSize = 0.05; 
const float _Density = 0.1;
const float _OpacityThreshold = 0.95;

const float _SphereRadius = 1.2;
const float _NoiseFreq = 0.5;
const float _NoiseAmp = 2.0;

const vec4 innerColor = vec4(0.7, 0.7, 0.7, _Density);
const vec4 outerColor = vec4(1.0, 1.0, 1.0, 0.0);

//const vec3 sunDir = vec3(-0.577, 0.577, 0.577);
const vec3 sunDir = vec3(-0.666, 0.333, 0.666);

// Description : Array and textureless GLSL 2D/3D/4D simplex 
//               noise functions.
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

vec4 mod289(vec4 x) {
  return x - floor(x * (1.0 / 289.0)) * 289.0;
}

vec4 permute(vec4 x) {
     return mod289(((x*34.0)+1.0)*x);
}

vec4 taylorInvSqrt(vec4 r)
{
  return 1.79284291400159 - 0.85373472095314 * r;
}

float snoise(vec3 v)
  { 
  const vec2  C = vec2(1.0/6.0, 1.0/3.0) ;
  const vec4  D = vec4(0.0, 0.5, 1.0, 2.0);

  // First corner
  vec3 i  = floor(v + dot(v, C.yyy) );
  vec3 x0 =   v - i + dot(i, C.xxx) ;

  // Other corners
  vec3 g = step(x0.yzx, x0.xyz);	  
  vec3 l = 1.0 - g;
  vec3 i1 = min( g.xyz, l.zxy );
  vec3 i2 = max( g.xyz, l.zxy );

  //   x0 = x0 - 0.0 + 0.0 * C.xxx;
  //   x1 = x0 - i1  + 1.0 * C.xxx;
  //   x2 = x0 - i2  + 2.0 * C.xxx;
  //   x3 = x0 - 1.0 + 3.0 * C.xxx;
  vec3 x1 = x0 - i1 + C.xxx;
  vec3 x2 = x0 - i2 + C.yyy; // 2.0*C.x = 1/3 = C.y
  vec3 x3 = x0 - D.yyy;      // -1.0+3.0*C.x = -0.5 = -D.y

  // Permutations
  i = mod289(i); 
  vec4 p = permute( permute( permute( 
             i.z + vec4(0.0, i1.z, i2.z, 1.0 ))
           + i.y + vec4(0.0, i1.y, i2.y, 1.0 )) 
           + i.x + vec4(0.0, i1.x, i2.x, 1.0 ));

  // Gradients: 7x7 points over a square, mapped onto an octahedron.
  // The ring size 17*17 = 289 is close to a multiple of 49 (49*6 = 294)
  float n_ = 0.142857142857; // 1.0/7.0
  vec3  ns = n_ * D.wyz - D.xzx;

  vec4 j = p - 49.0 * floor(p * ns.z * ns.z);  //  mod(p,7*7)

  vec4 x_ = floor(j * ns.z);
  vec4 y_ = floor(j - 7.0 * x_ );    // mod(j,N)

  vec4 x = x_ *ns.x + ns.yyyy;
  vec4 y = y_ *ns.x + ns.yyyy;
  vec4 h = 1.0 - abs(x) - abs(y);

  vec4 b0 = vec4( x.xy, y.xy );
  vec4 b1 = vec4( x.zw, y.zw );

  //vec4 s0 = vec4(lessThan(b0,0.0))*2.0 - 1.0;
  //vec4 s1 = vec4(lessThan(b1,0.0))*2.0 - 1.0;
  vec4 s0 = floor(b0)*2.0 + 1.0;
  vec4 s1 = floor(b1)*2.0 + 1.0;
  vec4 sh = -step(h, vec4(0.0));

  vec4 a0 = b0.xzyw + s0.xzyw*sh.xxyy ;
  vec4 a1 = b1.xzyw + s1.xzyw*sh.zzww ;

  vec3 p0 = vec3(a0.xy,h.x);
  vec3 p1 = vec3(a0.zw,h.y);
  vec3 p2 = vec3(a1.xy,h.z);
  vec3 p3 = vec3(a1.zw,h.w);

  //Normalise gradients
  vec4 norm = taylorInvSqrt(vec4(dot(p0,p0), dot(p1,p1), dot(p2, p2), dot(p3,p3)));
  p0 *= norm.x;
  p1 *= norm.y;
  p2 *= norm.z;
  p3 *= norm.w;

  // Mix final noise value
  vec4 m = max(0.6 - vec4(dot(x0,x0), dot(x1,x1), dot(x2,x2), dot(x3,x3)), 0.0);
  m = m * m;
  return 42.0 * dot( m*m, vec4( dot(p0,x0), dot(p1,x1), 
                                dot(p2,x2), dot(p3,x3) ) );
}


float fbm(vec3 p)
{
    float f;
    f = 0.5000*snoise( p ); p = p*2.02;
    f += 0.2500*snoise( p ); p = p*2.03;
    f += 0.1250*snoise( p ); p = p*2.01;
    f += 0.0625*snoise( p );
    return f;
}

float fbm2(vec3 p)
{
    const int octaves = 4;
    float amp = 0.5;
    float freq = 1.0;
    float n = 0.0;	
    for(int i=0; i<octaves; i++) {
        n += snoise(p*freq)*amp;
	freq *= 2.1;
	amp *= 0.5;
    }
    return n;
}

// returns signed distance to surface
float distanceFunc(vec3 p)
{	
	p.x -= iGlobalTime;		// translate with time
	//p += snoise(p*0.5)*1.0;	// domain warp!
	
	vec3 q = p;	
	// repeat on grid
	q.xz = mod(q.xz - vec2(2.5), 5.0) - vec2(2.5);
    q.y *= 2.0;	// squash in y
	float d = length(q) - _SphereRadius;	// distance to sphere

	// offset distance with noise
	//p = normalize(p) * _SphereRadius;	// project noise point to sphere surface
	p.y -= iGlobalTime*0.3;	// animate noise with time
	d += fbm(p*_NoiseFreq) * _NoiseAmp;
	return d;
}

// map distance to color
vec4 shade(float d)
{	
	return mix(innerColor, outerColor, smoothstep(0.5, 1.0, d));
}

// maps position to color
vec4 volumeFunc(vec3 p)
{
	float d = distanceFunc(p);
	vec4 c = shade(d);
	c.rgb *= smoothstep(-1.0, 0.0, p.y)*0.5+0.5;	// fake shadows
	float r = length(p)*0.04;
	c.a *= exp(-r*r);	// fog
	return c;
}

vec3 sky(vec3 v)
{
	// gradient
	vec3 c = mix(vec3(0.0, 0.5, 1.0), vec3(0, 0.25, 0.5), abs(v.y));
	//vec3 c = mix(vec3(1.0, 0.5, 0.0), vec3(0, 0.5, 1.0), abs(sqrt(v.y)));
	float sun = pow(dot(v, sunDir), 200.0);
	c += sun*vec3(3.0, 2.0, 1.0);
	return c;
}

float sampleLight(vec3 pos)
{
	const int _LightSteps = 8;
	const float _ShadowDensity = 1.0;
	vec3 lightStep = (sunDir * 2.0) / float(_LightSteps);
	float t = 1.0;	// transmittance
	for(int i=0; i<_LightSteps; i++) {
		vec4 col = volumeFunc(pos);
		t *= max(0.0, 1.0 - col.a * _ShadowDensity);
		//if (t < 0.01)
			//break;
		pos += lightStep;
	}
	return t;
}

// ray march volume
vec4 rayMarch(vec3 rayOrigin, vec3 rayStep, vec4 sum, out vec3 pos)
{
	pos = rayOrigin;
	for(int i=0; i<_VolumeSteps; i++) {
		vec4 col = volumeFunc(pos);
#if 0
		// volume shadows
		if (col.a > 0.0) {
			col.rgb *= sampleLight(pos);		
		}
#endif		
		
#if 0
		sum = mix(sum, col, col.a);	// under operator for back-to-front
#else	
		col.rgb *= col.a;		// pre-multiply alpha
		sum = sum + col*(1.0 - sum.a);	// over operator for front-to-back
#endif
		
#if 0
		// exit early if opaque
        	if (sum.a > _OpacityThreshold)
            		break;
#endif		
		pos += rayStep;
		//rayStep *= 1.01;
	}
	return sum;
}

bool
intersectBox(vec3 ro, vec3 rd, vec3 boxmin, vec3 boxmax, out float tnear, out float tfar)
{
	// compute intersection of ray with all six bbox planes
	vec3 invR = 1.0 / rd;
	vec3 tbot = invR * (boxmin - ro);
	vec3 ttop = invR * (boxmax - ro);
	// re-order intersections to find smallest and largest on each axis
	vec3 tmin = min (ttop, tbot);
	vec3 tmax = max (ttop, tbot);
	// find the largest tmin and the smallest tmax
	vec2 t0 = max (tmin.xx, tmin.yz);
	tnear = max (t0.x, t0.y);
	t0 = min (tmax.xx, tmax.yz);
	tfar = min (t0.x, t0.y);
	// check for hit
	bool hit;
	if ((tnear > tfar)) 
		hit = false;
	else
		hit = true;
	return hit;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 p = (fragCoord.xy / iResolution.xy)*2.0-1.0;
    p.x *= iResolution.x / iResolution.y;
	
    float rotx = 2.5 + (iMouse.y / iResolution.y)*4.0;
    float roty =  -0.2 - (iMouse.x / iResolution.x)*4.0;

    float zoom = 4.0;

    // camera
    vec3 ro = zoom*normalize(vec3(cos(roty), cos(rotx), sin(roty)));
	
    vec3 ww = normalize(vec3(0.0,0.0,0.0) - ro);
    vec3 uu = normalize(cross( vec3(0.0,1.0,0.0), ww ));
    vec3 vv = normalize(cross(ww, uu));
    vec3 rd = normalize( p.x*uu + p.y*vv + 1.5*ww );
	
    // box
    vec3 boxMin = vec3(-50.0, 2.0, -50);
    vec3 boxMax = vec3(50.0, -2.0, 50);
    //vec3 boxMin = vec3(-3.0, -2.0, -3.0);
    //vec3 boxMax = vec3(3.0, 2.0, 3.0);
	
    float tnear, tfar;
    bool hit = intersectBox(ro, rd, boxMin, boxMax, tnear, tfar);
    tnear = max(tnear, 0.0);
    tfar = max(tfar, 0.0);	

	vec3 pnear = ro+rd*tnear;
    vec3 pfar = ro+rd*tfar;
	
    //ro = pfar; rd = -rd; // back to front
    ro = pnear;	// front to back
    float stepSize = length(pfar - pnear) / float(_VolumeSteps);	
	
    vec3 hitPos;
    //vec4 col = vec4(0, 0.25, 0.5, 0);
    //vec4 col = vec4(sky(rd), 0);
    vec4 col = vec4(0);
    if (hit) {
	    //col = rayMarch(ro, rd*_StepSize, col, hitPos);
	    col = rayMarch(ro, rd*stepSize, col, hitPos);
    }

    // blend sun under clouds
    col += vec4(sky(rd), 0)*(1.0 - col.w);

    //col *= smoothstep(4.0, 0.7, dot(p, p));
	
    fragColor = vec4(col.rgb, 1.0);
}