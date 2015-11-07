// ----------------------------------------------------------------------------
// Main Rendering function
// ----------------------------------------------------------------------------

void mainImage(
	_out(vec4) fragColor,
	_in(vec2) fragCoord)
{
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
	//      R   ^ +Y                  ^ +Y             E eye/ray origin
	//       .  |\                    |     . R        R primary ray
	//         .| \                   |   .            @ fov angle
	//   -Z     | .\   +Z             | .
	//    ------0---E--->   +X -------0-------> -X
	//          | @/                  |
	//          | /                   |
	//          |/                    | -Y
	//           -Y
	//
	// NOTE: everything is expressed in this space, NOT world

	// assuming screen width is larger than height 
	vec2 aspect_ratio = vec2(u_res.x / u_res.y, 1);
	// field of view
	float fov = tan(radians(30.0));

	// antialising
#if 0
#define MSAA_PASSES 4
	float offset = 0.25;
	float ofst_x = offset * aspect_ratio.x;
	float ofst_y = offset;
	vec2 msaa[MSAA_PASSES];
	msaa[0] = vec2(-ofst_x, -ofst_y);
	msaa[1] = vec2(-ofst_x, +ofst_y);
	msaa[2] = vec2(+ofst_x, -ofst_y);
	msaa[3] = vec2(+ofst_x, +ofst_y);
#else
#define MSAA_PASSES 1
	vec2 msaa[MSAA_PASSES];
	msaa[0] = vec2(.5, .5);
#endif

	vec3 color = vec3(0, 0, 0);

	vec3 eye, look_at;
	setup_camera(eye, look_at);

	setup_scene();

	for (int i = 0; i < MSAA_PASSES; i++) {
		vec2 point_ndc = (fragCoord.xy + msaa[i]) / u_res.xy;
#ifdef HLSL
		point_ndc.y = 1. - point_ndc.y;
#endif
		vec3 point_cam = vec3(
			(2.0 * point_ndc - 1.0) * aspect_ratio * fov,
			-1.0);

		ray_t ray = get_primary_ray(point_cam, eye, look_at);

		color += render(ray) / float(MSAA_PASSES);
	}

	fragColor = vec4(corect_gamma(color, 2.25), 1);
}
