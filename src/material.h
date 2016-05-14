// ----------------------------------------------------------------------------
// Materials system definitions
// ----------------------------------------------------------------------------

struct material_t {
	vec3 base_color;
	float metallic;
	float roughness;
	float ior; // index of refraction
	float reflectivity;
	float translucency;
};

#define num_materials 8
#define mat_invalid -1
#define mat_debug 0
_mutable(material_t) materials[num_materials];

material_t get_material(
	_in(int) index
){
	material_t mat;

	for (int i = 0; i < num_materials; ++i) {
		if (i == index) {
			mat = materials[i];
			break;
		}
	}

	return mat;
}

