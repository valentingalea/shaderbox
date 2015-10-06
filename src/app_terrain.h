// Cloulds
// https://upload.wikimedia.org/wikipedia/commons/thumb/5/57/Cloud_types_en.svg/960px-Cloud_types_en.svg.png
// http://oceanservice.noaa.gov/education/yos/resource/JetStream/synoptic/clouds_max.htm#max
// http://www.metoffice.gov.uk/learning/clouds/cloud-spotting-guide
// http://www.srh.noaa.gov/srh/jetstream/clouds/cloudwise/types.html

// TODO: replace
float hash( float n ) { return fract(sin(n)*753.5453123); }
float noise( in vec3 x )
{
    vec3 p = floor(x);
    vec3 f = fract(x);
    f = f*f*(3.0-2.0*f);
	
    float n = p.x + p.y*157.0 + 113.0*p.z;
    return mix(mix(mix( hash(n+  0.0), hash(n+  1.0),f.x),
                   mix( hash(n+157.0), hash(n+158.0),f.x),f.y),
               mix(mix( hash(n+113.0), hash(n+114.0),f.x),
                   mix( hash(n+270.0), hash(n+271.0),f.x),f.y),f.z);
}
float fbm(vec3 pos)
{
    float f = 0.;
    vec3 q = pos;
    f  = 0.5000*noise( q ); q = q*2.01;
    f += 0.2500*noise( q ); q = q*2.02;
    f += 0.1250*noise( q ); q = q*2.03;
    f += 0.0625*noise( q ); q = q*2.01;    
    return f;
}

// GLSL/C++ compatibility layer
#ifdef __cplusplus
#define _in(T) const T &
#define _inout(T) T &
#define _out(T) T &
#define _begin {
#define _end }
#else
#define _in(T) const in T
#define _inout(T) inout T
#define _out(T) out T
#define _begin (
#define _end )
#endif

#define PI 3.14159265359

// Shadertoy specific uniforms
#define u_res iResolution
#define u_time iGlobalTime
#define u_mouse iMouse

struct ray_t {
	vec3 origin;
	vec3 direction;
};

ray_t get_primary_ray(_in(vec3) cam_local_point, _inout(vec3) cam_origin, _inout(vec3) cam_look_at)
{
	vec3 fwd = normalize(cam_look_at - cam_origin);
	vec3 up = vec3(0, 1, 0);
	vec3 right = cross(up, fwd);
	up = cross(fwd, right);

	return ray_t _begin
		cam_origin,
		normalize(fwd + up * cam_local_point.y + right * cam_local_point.x)
	_end;
}

mat3 rotate_around_x(_in(float) angle_degrees)
{
	float angle = radians(angle_degrees);
	float _sin = sin(angle);
	float _cos = cos(angle);
	return mat3(1, 0, 0, 0, _cos, -_sin, 0, _sin, _cos);
}

vec3 sun_dir = vec3(0, 1, 0);

vec3 render_sky_color(_in(ray_t) eye)
{
    float sun_dot = clamp(dot(eye.direction, sun_dir), 0., 1.);
    
    // colour scheme taken from https://www.shadertoy.com/view/MlSSR1
	vec3 blue = vec3(0.3,.55,0.8);
    vec3 red = vec3(0.8,0.8,0.6);        
    vec3 sky = mix(blue, red, 2.45*pow(sun_dot, 228.));
    
    return sky;// * (1. - 0.18*eye.direction);
}

vec3 render_clouds(_in(ray_t) eye)
{
	const int steps = 5;
	const float thickness = 20;
	
	float march_dist = 0.;
	float march_step = thickness / float (steps);
	
    // create vanishing point by projection with Y
    vec3 proj = eye.direction / (eye.direction.y + .225);
    
    // colors
    vec4 src = vec4 (0);
    vec4 dst = vec4 (0);
    
    for (int i = 0; i < steps; i++) {
    	vec3 sample = eye.origin + 
    	proj *
    	(march_dist + march_step * .5);
    	
    	// density func
    	// layer of cloulds + coverage
    	float n = fbm(sample * .242131 + u_time);
		src = vec4 (smoothstep(.33, 1., n));
		
		// decay func (could be beer law)
		src.a *= 0.5;

		// composition - front to back
		src.rgb *= src.a;
		dst = (1.0f - dst.a)*src + dst;
		
		// early out of opaque
		if (dst.a > .995) break;

    	march_dist += march_step;
    }
    
    // add horizon (hide lower artifact/reflection)
    // linear interp the Y with pow func that 
    // ramps up fast at the end, 0 otherwise
    dst = mix(dst, vec4 (0.), pow(1. - max(eye.direction.y, 0.), 8.));
    
    return dst.rgb;
}

void mainImage(_out(vec4) fragColor, _in(vec2) fragCoord)
{
	vec2 aspect_ratio = vec2(u_res.x / u_res.y, 1);
	float fov = tan(radians(45.0));
	vec2 point_ndc = fragCoord.xy / u_res.xy;
	vec3 point_cam = vec3((2.0 * point_ndc - 1.0) * aspect_ratio * fov, -1.0);

	vec3 col = vec3(0);

	//mat3 rot = rotate_around_x(-abs(sin(u_time / 2.)) * 90.);
	//sun_dir *= rot;
    
    vec3 eye = vec3 (0, 1., 0);
    vec3 look_at = vec3 (0, 1.5, -1);
	ray_t eye_ray = get_primary_ray(point_cam, eye, look_at);
    
    col += render_sky_color(eye_ray);
    col += render_clouds(eye_ray);
      
    fragColor = vec4(col, 1);
}