// ----------------------------------------------------------------------------
// Windows / DirectX specific
// ----------------------------------------------------------------------------
#include <cstdio>
#include <memory>
#include <new>
#include <thread>

typedef unsigned long DWORD;
#include "../../../lib/DirectXTex/DirectXTex/DDS.h"
using namespace DirectX;

struct DDS
{
	DWORD dwMagic;
	DDS_HEADER header;
	DDS_HEADER_DXT10 header10;
};

// ----------------------------------------------------------------------------
// CxxSwizzle support
// ----------------------------------------------------------------------------
#pragma warning(disable: 4244) // disable return implicit conversion warning
#pragma warning(disable: 4305) // disable truncation warning

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

// ----------------------------------------------------------------------------
// GLSL layer
// ----------------------------------------------------------------------------
#include "../../../src/def.h"

#include "../../../lib/ashima-noise/src/common.glsl"
#include "../../../lib/ashima-noise/src/classicnoise3d.glsl"
#include "../../../lib/ashima-noise/src/noise3d.glsl"
#include "../../../lib/ashima-noise/src/cellular3d.glsl"
#include "../../../src/noise_worley.h"

#include "../../../src/fbm.h"
//DECL_FBM_FUNC_TILE(fbm_dds, 4, (1. - noise_w(p, L).r))
DECL_FBM_FUNC_TILE(fbm_dds, 4, abs(pcnoise(p, L)))

// ----------------------------------------------------------------------------
// Main
// ----------------------------------------------------------------------------
int main(int argc, char* argv[])
{
	constexpr size_t size = 128;
	constexpr size_t channels = 1;
	using FLOAT = float;

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

	const size_t total_size = size * size* size * (sizeof(FLOAT) * channels);
	std::unique_ptr<FLOAT[]> data;
	data.reset(new (std::nothrow) FLOAT[total_size]);
	if (!data) {
		return 1;
	}

	auto worker = [&](size_t start, size_t count) {
		for (size_t z = start; z < start + count; z++) {
			for (size_t y = 0; y < size; y++) {
				for (size_t x = 0; x < size; x++) {
					vec3 input = (vec3(x, y, z) + .5f) / FLOAT(size);
					*(data.get() + size*size*z + size*y + x) =
						fbm_dds(input * 2.f, 2.f, .5f, .5);
				}
			}
		}
		printf("...finished [%i..%i]\n", start, start + count);
	};

	printf("starting work...\n");
	constexpr size_t size_quota = size / 4;
	std::thread w1{ worker, size_quota * 0, size_quota };
	std::thread w2{ worker, size_quota * 1, size_quota };
	std::thread w3{ worker, size_quota * 2, size_quota };
	std::thread w4{ worker, size_quota * 3, size_quota };
	w1.join();
	w2.join();
	w3.join();
	w4.join();

	auto file_closer = [](FILE* f) { fclose(f); };
	std::unique_ptr<FILE, decltype(file_closer)> file = { fopen("noise3d.dds", "wb+"), file_closer };
	if (!file) {
		return 2;
	}
	fwrite(&dds, sizeof(dds), 1, file.get());
	fwrite(data.get(), sizeof(FLOAT) * channels, size * size * size, file.get());

	return 0;
}