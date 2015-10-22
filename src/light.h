// ----------------------------------------------------------------------------
// Lighting models
// ----------------------------------------------------------------------------

#define LIGHT_POINT 1
#define LIGHT_DIR 2

struct light_t {
	int type;
	vec3 L; // origin for point lights, direction otherwise
	vec3 color;
};

light_t lights[8];

vec3 ambient_light = vec3(.01, .01, .01);

vec3 get_light_direction(
	_in(light_t) light,
	_in(hit_t) P
){
	if (light.type == LIGHT_DIR) {
		return light.L;
	} else {
		return normalize(light.L - P.origin);
	}
}

//     R       V    N    H      L         L dir to light       
//      ^      ^    ^    ^     ^          V dir to eye
//        .     \   |   /    .            N normal
//          .    \  |  /   .              H half between L and V
//            .   \ | /  .                R reflected
//  n1          .  \|/ .                  O hit point             
// -----------------O----------------     T refracted
//  n2             .                      n1 index of refraction of outgoing medium
//                .                       n2 index of refraction of incoming medium
//               .
//              .
//             .
//           \/_ T
//

vec3 illum_blinn_phong(
	_in(vec3) V,
	_in(vec3) L,
	_in(hit_t) hit,
	_in(material_t) mat
){
	vec3 diffuse = max(0., dot(L, hit.normal)) * mat.base_color;

	float spec_factor = 50.;
#if 0 // Blinn specular
	vec3 H = normalize(L + V);
	vec3 specular = pow(max(0., dot(H, hit.normal)), spec_factor); // * light.color * specular color
#else // Phong specular
	vec3 R = reflect(-L, hit.normal);
	vec3 specular = pow(max(0., dot(R, V)), spec_factor) * vec3(1, 1, 1); // * light.color * specular color
#endif

	return diffuse + specular;
}

vec3 illum_cook_torrance(
	_in(vec3) V,
	_in(vec3) L,
	_in(hit_t) hit,
	_in(material_t) mat
){
	vec3 H = normalize(L + V);
	float NdotL = dot(hit.normal, L);
	float NdotH = dot(hit.normal, H);
	float NdotV = dot(hit.normal, V);
	float VdotH = dot(V, H);

	// geometric term
	float geo_a = (2. * NdotH * NdotV) / VdotH;
	float geo_b = (2. * NdotH * NdotL) / VdotH;
	float geo_term = min(1., min(geo_a, geo_b));

	// roughness term -using Beckmann Distribution
	float rough_sq = mat.roughness * mat.roughness;
	float rough_a = 1. / (rough_sq * NdotH * NdotH * NdotH * NdotH);
	float rough_exp = (NdotH * NdotH - 1.) / (rough_sq * NdotH * NdotH);
	float rough_term = rough_a * exp(rough_exp);

	// Fresnel term
	float fresnel_term = fresnel_factor(1., mat.ior, VdotH);

	float specular = (geo_term * rough_term * fresnel_term) / (PI * NdotV * NdotL);
	return max(0., NdotL) * (specular + mat.base_color);
}
