// ----------------------------------------------------------------------------
// Various 3D utilities functions
// ----------------------------------------------------------------------------

ray_t get_primary_ray(
	_in(vec3) cam_local_point,
	_inout(vec3) cam_origin,
	_inout(vec3) cam_look_at
){
	vec3 fwd = normalize(cam_look_at - cam_origin);
	vec3 up = vec3(0, 1, 0);
	vec3 right = cross(up, fwd);
	up = cross(fwd, right);

	ray_t r = _begin(ray_t)
		cam_origin,
		normalize(fwd + up * cam_local_point.y + right * cam_local_point.x)
		_end;
	return r;
}

mat3 rotate_around_z(_in(float) angle_degrees)
{
	float angle = radians(angle_degrees);
	float _sin = sin(angle);
	float _cos = cos(angle);
	return mat3(_cos, -_sin, 0, _sin, _cos, 0, 0, 0, 1);
}

mat3 rotate_around_y(_in(float) angle_degrees)
{
	float angle = radians(angle_degrees);
	float _sin = sin(angle);
	float _cos = cos(angle);
	return mat3(_cos, 0, _sin, 0, 1, 0, -_sin, 0, _cos);
}

mat3 rotate_around_x(_in(float) angle_degrees)
{
	float angle = radians(angle_degrees);
	float _sin = sin(angle);
	float _cos = cos(angle);
	return mat3(1, 0, 0, 0, _cos, -_sin, 0, _sin, _cos);
}

vec3 corect_gamma(_in(vec3) color, _in(float) gamma)
{
	float p = 1.0 / gamma;
	return vec3(pow(color.r, p), pow(color.g, p), pow(color.b, p));
}

#ifdef __cplusplus
vec3 faceforward(_in(vec3) N, _in(vec3) I, _in(vec3) Nref)
{
	return dot(Nref, I) < 0 ? N : -N;
}
#endif
