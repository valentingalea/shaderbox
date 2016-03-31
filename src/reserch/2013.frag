//uniform float g_time;
varying vec2 p;

float t = iGlobalTime * .5;
const float c_pi = 3.1415926;

struct ray_t { vec3 origin, dir; };
#define RT(r,l) (r.origin + r.dir * l)
struct material_t { vec3 diffuse; float specular; };
struct light_dir_t { vec3 color, dir; };

material_t cmat = material_t(vec3(.6, .75, .2), 500000.);
light_dir_t sun = light_dir_t(vec3(1., .9, .5), normalize(vec3(1.3, 1., -.2)));

// Make a ray using normalized pixel position, eye position and focus point
ray_t lookAtDir(in vec3 uv_dir, in vec3 pos, in vec3 at) {
	vec3 f = normalize(at - pos);
	vec3 r = cross(f, vec3(0.,1.,0.));
	vec3 u = cross(r, f);
	return ray_t(pos, normalize(uv_dir.x * r + uv_dir.y * u + uv_dir.z * f));
}

float hash(in float x) { return fract(sin(x*.0007)*29835.24389); }
float hash(in vec2 x) { return hash(dot(x,vec2(23.17,17.23))); }
float hash(in vec3 x) { return hash(dot(x,vec3(23.17,17.23,42.5))); }
float noise_value(in vec2 p) {
	vec2 e=vec2(1.,0.), F = floor(p), f = fract(p), k = (3. - 2.*f) * f * f;
	return mix(mix(hash(F),      hash(F+e.xy), k.x),
			   mix(hash(F+e.yx), hash(F+e.xx), k.x), k.y);
}
float noise_value(in vec3 p) {
	vec2 e = vec2(1.,0.);
	vec3 F = floor(p), f = fract(p), k = (3. - 2.*f) * f * f;
	return mix(
		   mix(mix(hash(F),       hash(F+e.xyy), k.x),
			   mix(hash(F+e.yxy), hash(F+e.xxy), k.x), k.y),
		   mix(mix(hash(F+e.yyx), hash(F+e.xyx), k.x),
			   mix(hash(F+e.yxx), hash(F+e.xxx), k.x), k.y), k.z);
}

float terrain_height(in vec3 at) {
	float k = .5;
	vec2 p = at.xz * .0005;
	mat2 m = mat2(1.3,  1.5, 1.5, -1.3);
	float h = 0.;//noise_value(p) * k;
	for (int i = 0; i < 8; ++i) { h += noise_value(p) * k; k *= .5; p = m * p; }
	return 1200. * h * h * h * h;
}
	
float terrain_distance(in vec3 at) {
	return min(1.+at.y,at.y - terrain_height(at));
}

float cloud_noise(in vec3 p) {
	float k = .5, v = 0.;
    v += noise_value(p) * k; k *= .5; p =  p * 2.01;
    v += noise_value(p) * k; k *= .5; p =  p * 2.01;
    v += noise_value(p) * k; k *= .5; p =  p * 2.01;
    v += noise_value(p) * k; k *= .5; p =  p * 2.01;
	return v;
}

#define MAX_DISTANCE 23000.
#define CLOUDS_HEIGHT_MAX 3000.
#define CLOUDS_HEIGHT_MIN 1000.

vec4 accum = vec4(0.);
float prevV = 0.;
float dist = 0.;
float solid = 0.;

float world(in vec3 at, in vec3 dir) {
	if (at.y > CLOUDS_HEIGHT_MAX)
		return (dir.y < 0.) ? (at.y - CLOUDS_HEIGHT_MAX + 10.) / -dir.y : 1000.;
	if (at.y > CLOUDS_HEIGHT_MIN) {
		vec3 cs = at*.001 - t *.3;
		float k =0.,// sin((at.y - CLOUDS_HEIGHT_MIN) / (CLOUDS_HEIGHT_MAX - CLOUDS_HEIGHT_MIN) * c_pi),
			v = cloud_noise(cs) * k,
			V = max(0., v * v - .13);
		if (V > .0) {
			float A = ((V + prevV) * .5 * dist) * .005;
			//vec3 cloudcolor = vec3(.125);// * V;// + sun.color * max(0.,(-cloud_noise(cs+sun.dir*.01)+v))*7.725;
			vec3 cloudcolor = vec3(1.) + sun.color * max(0.,(-cloud_noise(cs+sun.dir*100.)+v)) * 2.;
			accum.rgb += cloudcolor * A * (1. - accum.a);
			accum.a += A;
		}
		prevV = V;
		return 1. + 200. * (1. - V);
	}
	float H = terrain_distance(at);
	return (dir.y <= 0.) ? H : min(H, (CLOUDS_HEIGHT_MIN + 10. - at.y) / dir.y);
}

vec3 normal(in vec3 at) {
	vec2 e = vec2(.1, .0);
	return normalize(vec3(world(at+e.xyy, vec3(0.))-world(at, vec3(0.)),
		  world(at+e.yxy, vec3(0.))-world(at, vec3(0.)),
		  world(at+e.yyx, vec3(0.))-world(at, vec3(0.))));
}

float trace(in ray_t ray, in float maxPath) {
	float path = 0.;
	for (int i = 0; i < 128; ++i) {
		dist = world(RT(ray, path), ray.dir);
		if (dist < .001*path) { solid = 1.; break; }
		if (accum.a >= 1.) break;
		path += dist;
		if (path > maxPath) break;
	}
	return path;
}

vec3 background(in vec3 dir) {
    vec3 ac = vec3(101.,133.,162.) / 255.;
	return ac + sun.color * pow(max(0., dot(dir, sun.dir)), 40.);
}

vec3 light(in vec3 at, in vec3 from, in light_dir_t l) {
	vec3 n = normal(at), h = normalize(from + l.dir);
	//if (trace(ray_t(at+n*1., l.dir), 3000.) < 3000.) return vec3(0.);
	return l.color * (
		cmat.diffuse * max(0., dot(n, l.dir))
		+ pow(max(0., dot(n, h)), cmat.specular) * (cmat.specular + 8.) / (8. * c_pi));
}

void update_mat(in vec3 at) {
	float h = at.y;// + 30. * (noise_value(at.zx*.00002) - .5);
	if (h < 125.) {
		cmat.diffuse = vec3(.2, .4, .9);
		cmat.specular = 30.;
	} else if(h < 150.) {
		cmat.diffuse = vec3(.8, .8, .2);
	} else if(h < 220.) {
		cmat.diffuse = vec3(.8, .9, .2);
	} else if(h > 400.) {
		cmat.diffuse = vec3(1.);
	} else {
		cmat.specular = 100000.;
	}
}

void mainImage( out vec4 fragColor, in vec2 fragCoord ) {
	// Calculate normalized and aspect-corrected pixel position 
	float aspect = iResolution.x / iResolution.y;
	vec2 uv = (fragCoord.xy / iResolution.xy * 2. - 1.) * vec2(aspect, 1.);
	
	ray_t ray;

		vec2 mp = iMouse.xy / iResolution.xy;
		float a = 2. * c_pi * mp.x;
		float s = sin(a), c = cos(a);
		mat3 r = mat3(c, 0., -s, 0., 1., 0., s, 0., c);
		ray = lookAtDir(normalize(vec3(uv, 2.)),
		 			    vec3(1500.*t, 0. ,0.)+r*vec3(6000.,100.+6000.*(1.-mp.y),6000.),
						vec3(1500.*t, 0. ,0.));
	

	ray.origin.y = max(terrain_height(ray.origin)+10., ray.origin.y);
	float path = trace(ray, MAX_DISTANCE);

	vec3 color = vec3(0.);
	accum.a = min(accum.a, 1.);
// terrain coloring
	if (solid > 0.) {
		vec3 at = RT(ray,path);
		update_mat(at);
		color = light(at, -ray.dir, sun);
	} else {
		color = accum.rgb;
	}

// alpha blending of clouds
    float klen = clamp(pow(path / MAX_DISTANCE, 2.), 0., 1.);
	color = mix(color, background(ray.dir), klen) * (1.-accum.a) + accum.rgb;
    
// sun coloring
    //color += .64 * sun.color * pow(max(0., dot(ray.dir, sun.dir)), 2.);

	fragColor = vec4(color, 0.);
}