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
// VML support
// ----------------------------------------------------------------------------
#pragma warning(disable: 4244) // disable return implicit conversion warning
#pragma warning(disable: 4305) // disable truncation warning

#include <vml/vector.h>
#include <vml/matrix.h>
#include <vml/vector_functions.h>

using  vec4 = vml::vector<float, 0, 1, 2, 3>;
using  vec3 = vml::vector<float, 0, 1, 2>;
using  vec2 = vml::vector<float, 0, 1>;
using   _01 = vml::indices_pack<0, 1>;
using  _012 = vml::indices_pack<0, 1, 2>;
using _0123 = vml::indices_pack<0, 1, 2, 3>;
using  mat2 = vml::matrix<float, vml::vector, _01, _01>;
using  mat3 = vml::matrix<float, vml::vector, _012, _012>;
using  mat4 = vml::matrix<float, vml::vector, _0123, _0123>;

// ----------------------------------------------------------------------------
// GLSL layer
// ----------------------------------------------------------------------------
#include "../../../src/def.h"
#include "../../../src/util.h"
//#include "../../../lib/ashima-noise/src/common.glsl"
//#include "../../../lib/ashima-noise/src/classicnoise3d.glsl"
//#include "../../../lib/ashima-noise/src/noise3d.glsl"
//#include "../../../lib/ashima-noise/src/cellular3d.glsl"
#include "../../../src/noise_worley.h"
#include "../../../src/fbm.h"

DECL_FBM_FUNC_TILE(fbm_worley_tile, 4, (1. - (noise_w(p, L).r + .25)))
//DECL_FBM_FUNC_TILE(fbm_perlin_tile, 4, abs(pcnoise(p, L)))

float fbm_dds(vec3 &pos)
{
//	float p = fbm_perlin_tile(pos, 2., 1., .5);
//	float w = 1. - fbm_worley_tile(pos, 4., 1., .5);
//	return remap(p, -w, 1., 0., 1.);
	return fbm_worley_tile(pos, 2., 1., .5);
}

// ----------------------------------------------------------------------------
// Main
// ----------------------------------------------------------------------------
int main(int argc, char* argv[])
{
	constexpr size_t size = 128;
	constexpr size_t channels = 4;
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

	DXGI_FORMAT fmt[] = { DXGI_FORMAT_R32_FLOAT, DXGI_FORMAT_R32G32_FLOAT, DXGI_FORMAT_R32G32B32_UINT, DXGI_FORMAT_R32G32B32A32_FLOAT };
	dds.header10.dxgiFormat = fmt[channels - 1];
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
				auto ptr = data.get() + size*size*channels*z + size*channels*y;

				for (size_t x = 0; x < size; x++) {
					vec3 pos = (vec3(x, y, z) + .5f) / FLOAT(size);
					*ptr++ = fbm_dds(pos);
					*ptr++ = 0.f;
					*ptr++ = 0.f;
					*ptr++ = 0.f;
				}
			}
		}
		printf("...finished [%i..%i]\n", start, start + count);
	};

	printf("starting work...\n");
	using tclock_t = std::chrono::steady_clock;
	using tpoint_t = tclock_t::time_point;
	tpoint_t t_start = tclock_t::now();

	constexpr size_t size_quota = size / 4;
	std::thread w1{ worker, size_quota * 0, size_quota };
	std::thread w2{ worker, size_quota * 1, size_quota };
	std::thread w3{ worker, size_quota * 2, size_quota };
	std::thread w4{ worker, size_quota * 3, size_quota };
	w1.join();
	w2.join();
	w3.join();
	w4.join();

	tpoint_t t_end = tclock_t::now();
	using tdiff_t = std::chrono::duration<float, std::chrono::seconds::period>;
	auto elapsted = tdiff_t(t_end - t_start).count();
	printf("total running time: %f\n", elapsted);

#if 1
	time_t rawtime;
	time(&rawtime);
	tm timeinfo;
	localtime_s(&timeinfo, &rawtime);
	constexpr int SZ = 128;
	char dds_file[SZ];
	strftime(dds_file, SZ, "noise3d-%Y-%m-%d-%H-%M-%S.dds\0", &timeinfo);

	auto file_closer = [](FILE* f) { fclose(f); };
	std::unique_ptr<FILE, decltype(file_closer)> file = { fopen(dds_file, "wb+"), file_closer };
	if (!file) {
		return 2;
	}
	fwrite(&dds, sizeof(dds), 1, file.get());
	fwrite(data.get(), sizeof(FLOAT) * channels, size * size * size, file.get());
#endif

	return 0;
}