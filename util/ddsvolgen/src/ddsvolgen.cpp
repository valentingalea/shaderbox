#include <swizzle/glsl/naive/vector.h>
#include <swizzle/glsl/naive/matrix.h>
#include <swizzle/glsl/texture_functions.h>
#include <swizzle/glsl/vector_functions.h>
typedef float real_t;
typedef swizzle::glsl::naive::vector< int, 2 > ivec2;
typedef swizzle::glsl::naive::vector< real_t, 2 > vec2;
typedef swizzle::glsl::naive::vector< real_t, 3 > vec3;
typedef swizzle::glsl::naive::vector< real_t, 4 > vec4;
typedef swizzle::glsl::naive::matrix< swizzle::glsl::naive::vector, real_t, 2, 2> mat2;
typedef swizzle::glsl::naive::matrix< swizzle::glsl::naive::vector, real_t, 3, 3> mat3;
typedef swizzle::glsl::naive::matrix< swizzle::glsl::naive::vector, real_t, 4, 4> mat4;

#include "../../../src/def.h"

//#include "../../../src/noise_iq.h"
//#define noise(x) noise_iq(x)

//#include "../../../src/noise_worley.h"
//#define noise(x) (1. - noise_w(x).r)
//#define noise(x) abs( noise_iq(x / 8.) - (1. - (noise_w(x * 2.).r)))

#include "../../../lib/ashima-noise/src/classicnoise3D.glsl"
#define noise(x) cnoise(x)

#include "../../../src/fbm.h"
DECL_FBM_FUNC(fbm_dds, 5, abs(pnoise(p, vec3(lacunarity))))

#include <d3d11.h>
#include "../../../lib/DirectXTex/DirectXTex/DDS.h"
#include <cstdio>
#include <memory>
#include <new>

using namespace DirectX;

struct DDS
{
	DWORD dwMagic;
	DDS_HEADER header;
	DDS_HEADER_DXT10 header10;
};

int main(int argc, char* argv[])
{
	auto size = 128;
	auto channels = 1;

	DDS dds = { 0 };

	dds.dwMagic = DDS_MAGIC;

	dds.header.dwSize = sizeof(DDS_HEADER);
	dds.header.dwFlags = DDS_HEADER_FLAGS_TEXTURE | DDS_HEADER_FLAGS_VOLUME | DDS_HEADER_FLAGS_PITCH;
	dds.header.dwHeight = size;
	dds.header.dwWidth = size;
	dds.header.dwDepth = size;
	dds.header.dwPitchOrLinearSize = (size * (sizeof(FLOAT) * channels) + 7) / 8;
	dds.header.dwMipMapCount = 0;
	dds.header.ddspf = DDSPF_DX10;
	dds.header.dwCaps = DDS_SURFACE_FLAGS_TEXTURE | DDS_SURFACE_FLAGS_CUBEMAP;
	dds.header.dwCaps2 = DDS_FLAGS_VOLUME;

	dds.header10.dxgiFormat = DXGI_FORMAT_R32_FLOAT;
	dds.header10.resourceDimension = DDS_DIMENSION_TEXTURE3D;
	dds.header10.arraySize = 1;
	dds.header10.miscFlag = 0;
	dds.header10.miscFlags2 = 0;

	const auto total_size = size * size* size * (sizeof(FLOAT) * channels);
	std::unique_ptr<FLOAT[]> data;
	data.reset(new (std::nothrow) FLOAT[total_size]);
	if (!data) {
		return 1;
	}

	#pragma loop( hint_parallel(8) )
	for (int z = 0; z < size; z++) {
		#pragma loop( hint_parallel(8) )
		for (int y = 0; y < size; y++) {
			#pragma loop( hint_parallel(8) )
			for (int x = 0; x < size; x++) {
				vec3 input = vec3(x, y, z) / FLOAT(size);
				*(data.get() + size*size*z + size*y + x) = 
					fbm_dds(input * 2.03, 2.64, .5, .5);
			}
		}
	}

	auto file_closer = [](FILE* f) { fclose(f); };
	std::unique_ptr<FILE, decltype(file_closer)> file = { fopen("noise3d.dds", "wb+"), file_closer };
	if (!file) {
		return 2;
	}
	fwrite(&dds, sizeof(dds), 1, file.get());
	fwrite(data.get(), sizeof(FLOAT) * channels, size * size * size, file.get());

	return 0;
}