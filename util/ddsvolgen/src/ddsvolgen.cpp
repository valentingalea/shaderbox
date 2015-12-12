#include <d3d11.h>
#include "../../../lib/DirectXTex/DirectXTex/DDS.h"
#include <cstdio>
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

	dds.header10.dxgiFormat = DXGI_FORMAT_R32G32B32A32_FLOAT;
	dds.header10.resourceDimension = DDS_DIMENSION_TEXTURE3D;
	dds.header10.arraySize = 1;
	dds.header10.miscFlag = 0;
	dds.header10.miscFlags2 = 0;

	const auto total_size = size * size* size * (sizeof(FLOAT) * channels);
	auto *data = new (std::nothrow) FLOAT[total_size];
	if (!data) {
		return 1;
	}

	for (int z = 0; z < size; z++) {
		for (int y = 0; y < size; y++) {
			for (int x = 0; x < size; x++) {
				*(data + size*size*z + size*y + x) = 1.;
			}
		}
	}

	FILE *file = fopen("noise3d.dds", "wb+");
	if (!file) {
		return 2;
	}
	fwrite(&dds, sizeof(dds), 1, file);
	fwrite(data, sizeof(FLOAT) * channels, size * size * size, file);
	fclose(file);

	delete[] data;

	return 0;
}