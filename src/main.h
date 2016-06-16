// ----------------------------------------------------------------------------
// Main Rendering function
// depends on external defines: FOV
// ----------------------------------------------------------------------------

void mainImage(
	_out(vec4) fragColor,
#ifdef SHADERTOY
	vec2 fragCoord
#else
	_in(vec2) fragCoord
#endif
){
	// The pipeline transform
	//
	// 1. fragCoord is in raster space [0..resolution]
	// 2. convert to NDC [0..1] by dividing to the resolution
	// 3. convert to camera space:
	//  a. xy gets [-1, +1] by 2 * NDC - 1; z fixed at -1
	//  c. apply aspect & fov
	//  d. apply the look-at algoritm which will
	//     produce the 3 camera axis:
	//
	//      R   ^ -Y                  ^ -Y             E eye/ray origin
	//       .  |\                    |     . R        R primary ray
	//         .| \                   |   .            @ fov angle
	//   -Z     | .\   +Z             | .
	//    ------0---E--->   -X -------0-------> +X
	//          | @/                  |
	//          | /                   |
	//          |/                    | +Y
	//           +Y
	//
	// NOTE: everything is expressed in this space, NOT world

	// assuming screen width is larger than height 
	vec2 aspect_ratio = vec2(u_res.x / u_res.y, 1);

	vec3 eye, look_at;
	setup_camera(eye, look_at);

	setup_scene();

	vec2 point_ndc = fragCoord.xy / u_res.xy;
#ifdef HLSL
		point_ndc.y = 1. - point_ndc.y;
#endif
	vec3 point_cam = vec3(
		(2.0 * point_ndc - 1.0) * aspect_ratio * FOV,
		-1.0);

	ray_t ray = get_primary_ray(point_cam, eye, look_at);

	vec3 color = render(ray, point_cam);

	fragColor = vec4(linear_to_srgb(color), 1);
}