// from http://the-witness.net/news/2012/11/scopeexit-in-c11/
template <typename F>
struct ScopeExit {
	ScopeExit(F f) : f(f) {}
	~ScopeExit() { f(); }
	F f;
};
template <typename F>
ScopeExit<F> MakeScopeExit(F f) {
	return ScopeExit<F>(f);
};
#define STRING_JOIN2(arg1, arg2) DO_STRING_JOIN2(arg1, arg2)
#define DO_STRING_JOIN2(arg1, arg2) arg1 ## arg2
#define SCOPE_EXIT(code) \
    auto STRING_JOIN2(scope_exit_, __LINE__) = MakeScopeExit([&](){code;})

#define SafeRelease(T) if (T) { T->Release(); T = NULL; }

#include <windowsx.h>
#include <d3d11.h>
#include <D3Dcompiler.h>
#include <DirectXMath.h>
#include <CRTDBG.H>
#include <cstdio>
#include <chrono>
#include "../../../lib/DirectXTex/DirectXTex/DirectXTex.h"

#define IMGUI_DISABLE_OBSOLETE_FUNCTIONS
#include <imgui.h>
#include <examples\directx11_example\imgui_impl_dx11.h>
IMGUI_API LRESULT   ImGui_ImplDX11_WndProcHandler(HWND hWnd, UINT msg, WPARAM wParam, LPARAM lParam);

void ShowError(LPCSTR szErrMsg, ID3D10Blob* pExtraErrorMsg = NULL)
{
	constexpr int MAX = 4096;
	char message[MAX];
	sprintf_s(message, MAX, szErrMsg);
	if (pExtraErrorMsg != NULL)
		sprintf_s(message, MAX, "%s\n%s", szErrMsg, (LPCTSTR)pExtraErrorMsg->GetBufferPointer());
	MessageBox(0, message, "Error", MB_ICONERROR);
}

LPCTSTR lpszClassName = "tinyDX11";
LPCTSTR lpszAppName = "hlsltoy";
constexpr int WIDTH = 1280;
constexpr int HEIGHT = 720;

DXGI_SWAP_CHAIN_DESC SwapChainDesc =
{
	{ WIDTH, HEIGHT, { 0, 0 }, DXGI_FORMAT_R8G8B8A8_UNORM, DXGI_MODE_SCANLINE_ORDER_UNSPECIFIED, DXGI_MODE_SCALING_UNSPECIFIED, },  // w, h, refreshrate, format, scanline ordering, scaling
	{ 1, 0 }, // number of multisamples per pixel, quality level
	DXGI_USAGE_RENDER_TARGET_OUTPUT, 1, // buffer usage, buffer count
	0, // output window (must be supplied later)
	1, DXGI_SWAP_EFFECT_DISCARD, 0 // windowed, swap effect, flags
};

// from http://altdevblog.com/2011/08/08/an-interesting-vertex-shader-trick/
CHAR szVertexShader[] =
"float4 main(uint id : SV_VertexID) : SV_Position {"
	"float2 tex = float2((id << 1) & 2, id & 2);"
	"return float4(tex * float2(2, -2) + float2(-1, 1), 0, 1);"
"}";

// from https://msdn.microsoft.com/en-us/library/windows/desktop/ff476521%28v=vs.85%29.aspx
ID3D11Texture2D* CreateTextureCheckboard(ID3D11Device *pd3dDevice, UINT w, UINT h, UINT checkFreq)
{
	UINT *buffer = new UINT[w * h];
	_ASSERT(buffer);
	SCOPE_EXIT(delete[] buffer);

	for (UINT y = 0; y < h; y++)
		for (UINT x = 0; x < w; x++)
			buffer[y*h + x] = ((x & checkFreq) == (y & checkFreq)) ? 0xff000000 : 0xffffffff;

	ID3D11Texture2D *tex = NULL;
	D3D11_TEXTURE2D_DESC tdesc = { w, h, 1 /*multisampled*/, 1 /*mip levels*/, DXGI_FORMAT_R8G8B8A8_UNORM, { 1, 0 } /*quality: no AA*/,
		D3D11_USAGE_DEFAULT, D3D11_BIND_SHADER_RESOURCE, 0 /*no CPU access*/, 0 /*flags*/ };
	D3D11_SUBRESOURCE_DATA tbsd = { buffer, w * 4 /*pitch*/, w * h * 4 /*pitch elem in array, unused*/ };

	HRESULT hr = pd3dDevice->CreateTexture2D(&tdesc, &tbsd, &tex);
	_ASSERT(SUCCEEDED(hr));

	return tex;
}

ID3D11ShaderResourceView* CreateNoiseTexture(ID3D11Device *pd3dDevice, LPCWSTR lpszPath)
{
	using namespace DirectX;

	TexMetadata mdata;
	HRESULT hr = GetMetadataFromDDSFile(lpszPath, DDS_FLAGS_NONE, mdata);
	_ASSERT(SUCCEEDED(hr));
	_ASSERT(mdata.IsVolumemap());

	ScratchImage image;
	hr = LoadFromDDSFile(lpszPath, DDS_FLAGS_NONE, &mdata, image);
	if (FAILED(hr)) {
		WCHAR buff[2048];
		swprintf_s(buff, L"Failed to load texture file\n\nFilename = %ls\nHRESULT %08X", lpszPath, hr);
		MessageBoxW(0, buff, L"Error", MB_ICONERROR);
		return NULL;
	}

	ID3D11ShaderResourceView* pView = NULL;
	hr = CreateShaderResourceView(pd3dDevice, image.GetImages(), image.GetImageCount(), mdata, &pView);
	_ASSERT(SUCCEEDED(hr));

	return pView;
}

POINT sMouseButtons = { 0, 0 };
POINT sMousePos = { 0, 0 };

LRESULT CALLBACK WndProc(HWND hWnd, UINT uMsg, WPARAM wParam, LPARAM lParam)
{
	if (uMsg == WM_CLOSE)
	{
		PostQuitMessage(0);
		return TRUE;
	}

	ImGui_ImplDX11_WndProcHandler(hWnd, uMsg, wParam, lParam);

	if (!ImGui::IsMouseHoveringAnyWindow())
	{
		switch (uMsg)
		{
		case WM_LBUTTONDOWN:
			sMouseButtons.x = TRUE;
			break;
		case WM_LBUTTONUP:
			sMouseButtons.x = FALSE;
			break;
		case WM_MOUSEMOVE:
			sMousePos.x = GET_X_LPARAM(lParam);
			sMousePos.y = GET_Y_LPARAM(lParam);
			break;
		}
	}

	return DefWindowProc(hWnd, uMsg, wParam, lParam);
}

int __stdcall WinMain(HINSTANCE hThisInstance, HINSTANCE hPrevInstance, LPSTR lpszArgument, int nCmdShow)
{
//
// Command line
//
	LPWSTR *szArglist = NULL;
	SCOPE_EXIT(if (szArglist) LocalFree(szArglist));
	int nArgs = 0;
	szArglist = CommandLineToArgvW(GetCommandLineW(), &nArgs);
	if (!szArglist || nArgs < 2)
	{
		ShowError("Minimal DX11 Shader Framework.\n\Usage:\nhlsltoy <shader file> [noise dds file] [noise dds file]");
		return ERROR_INVALID_COMMAND_LINE;
	}

//
// Win32 window
//
	HRESULT hr = NULL;
	HINSTANCE hinst = GetModuleHandle(NULL);
	WNDCLASS wc = { CS_HREDRAW | CS_VREDRAW | CS_OWNDC, (WNDPROC)WndProc, 0, 0, hinst, LoadIcon(NULL, IDI_WINLOGO), LoadCursor(NULL, IDC_ARROW), 0, 0, lpszClassName };
	RegisterClass(&wc);
	SCOPE_EXIT(UnregisterClass(lpszClassName, hinst));

	DWORD dwStyle = WS_POPUPWINDOW | WS_CAPTION | WS_MINIMIZEBOX | WS_VISIBLE;
	int sysWidth = ::GetSystemMetrics(SM_CXFULLSCREEN);
	int sysHeight = ::GetSystemMetrics(SM_CYFULLSCREEN);
	RECT desiredBox;
	SetRect(&desiredBox, 0, 0, WIDTH, HEIGHT);
	AdjustWindowRect(&desiredBox, dwStyle, FALSE);
	int newWidth = desiredBox.right - desiredBox.left;
	int newHeight = desiredBox.bottom - desiredBox.top;
	HWND hWnd = CreateWindow(lpszClassName, lpszAppName, dwStyle,
		sysWidth / 2 - newWidth / 2, sysHeight / 2 - newHeight / 2, newWidth, newHeight,
		0, 0, hinst, 0);
	SwapChainDesc.OutputWindow = hWnd;

//
// DX setting up device
//
	ID3D11Device* pd3dDevice = NULL;
	ID3D11DeviceContext* pImmediateContext = NULL;
	IDXGISwapChain* pSwapChain = NULL;
	SCOPE_EXIT(SafeRelease(pSwapChain));
	SCOPE_EXIT(SafeRelease(pImmediateContext));
	SCOPE_EXIT(SafeRelease(pd3dDevice));
	hr = D3D11CreateDeviceAndSwapChain(NULL, D3D_DRIVER_TYPE_HARDWARE, NULL,
#ifdef _DEBUG
		D3D11_CREATE_DEVICE_DEBUG,
#else
		D3D11_CREATE_DEVICE_SINGLETHREADED,
#endif
		0, 0, D3D11_SDK_VERSION, &SwapChainDesc, &pSwapChain, &pd3dDevice, NULL, &pImmediateContext);
	if (FAILED(hr)) return hr;

//
// Checkboard texture
//
	ID3D11Texture2D* pTex = CreateTextureCheckboard(pd3dDevice, 128, 128, 16);
	ID3D11ShaderResourceView *pTexV = NULL;
	SCOPE_EXIT(SafeRelease(pTex));
	SCOPE_EXIT(SafeRelease(pTexV));
	hr = pd3dDevice->CreateShaderResourceView(pTex, NULL/*whole res*/, &pTexV);
	_ASSERT(SUCCEEDED(hr));

//
// Noise textures
//
	ID3D11ShaderResourceView *pNoiseTexV = NULL;
	if (nArgs >= 3)
	{
		pNoiseTexV = CreateNoiseTexture(pd3dDevice, szArglist[2]);
	}
	SCOPE_EXIT(SafeRelease(pNoiseTexV));

	ID3D11ShaderResourceView *pNoiseTexV_2 = NULL;
	if (nArgs >= 4)
	{
		pNoiseTexV_2 = CreateNoiseTexture(pd3dDevice, szArglist[3]);
	}
	SCOPE_EXIT(SafeRelease(pNoiseTexV_2));

//
// Textures sampler
//
	ID3D11SamplerState *pSampler = NULL;
	SCOPE_EXIT(SafeRelease(pSampler));
	D3D11_SAMPLER_DESC sSamplerDesc = { D3D11_FILTER_MIN_MAG_MIP_LINEAR, D3D11_TEXTURE_ADDRESS_WRAP, D3D11_TEXTURE_ADDRESS_WRAP , D3D11_TEXTURE_ADDRESS_WRAP,
		0., 1, D3D11_COMPARISON_NEVER , {1, 1, 1, 1}, -FLT_MAX, FLT_MAX };
	hr = pd3dDevice->CreateSamplerState(&sSamplerDesc, &pSampler);
	_ASSERT(SUCCEEDED(hr));

//
// Backbuffer render target
//
	ID3D11Texture2D* pBackBuffer = NULL;
	SCOPE_EXIT(SafeRelease(pBackBuffer));
	hr = pSwapChain->GetBuffer(0, __uuidof(ID3D11Texture2D), (LPVOID*)&pBackBuffer);
	_ASSERT(SUCCEEDED(hr));

	ID3D11RenderTargetView* pRenderTargetView = NULL;
	SCOPE_EXIT(SafeRelease(pRenderTargetView));
	hr = pd3dDevice->CreateRenderTargetView(pBackBuffer, NULL, &pRenderTargetView);
	_ASSERT(SUCCEEDED(hr));

	pImmediateContext->OMSetRenderTargets(1, &pRenderTargetView, NULL);

//
// Wiewport & Rasterizer
//
	D3D11_VIEWPORT vp = { 0, 0, WIDTH, HEIGHT, 0., 1. };
	pImmediateContext->RSSetViewports(1, &vp);

	D3D11_RASTERIZER_DESC rasterizerDesc = { D3D11_FILL_SOLID, D3D11_CULL_NONE, FALSE, 0, 0., 0., FALSE, FALSE, FALSE, FALSE };
	ID3D11RasterizerState* pd3dRasterizerState = NULL;
	SCOPE_EXIT(SafeRelease(pd3dRasterizerState));
	hr = pd3dDevice->CreateRasterizerState(&rasterizerDesc, &pd3dRasterizerState);
	_ASSERT(SUCCEEDED(hr));

//
// Common tracking vars & settings
//
	ID3D10Blob* pErrorBlob = NULL;
	ID3D10Blob* pBlob = NULL;
	SCOPE_EXIT(SafeRelease(pErrorBlob));
	SCOPE_EXIT(SafeRelease(pBlob));
	UINT flags = D3DCOMPILE_ENABLE_STRICTNESS | D3DCOMPILE_IEEE_STRICTNESS | D3DCOMPILE_WARNINGS_ARE_ERRORS |
#ifdef _DEBUG
		D3DCOMPILE_DEBUG | D3DCOMPILE_SKIP_OPTIMIZATION
#else
		D3DCOMPILE_OPTIMIZATION_LEVEL3
#endif
		;

//
// Vertex shader
//
	ID3D11VertexShader* pVS = NULL;
	SCOPE_EXIT(SafeRelease(pVS));
	hr = D3DCompile(szVertexShader, sizeof(szVertexShader), 0, 0, 0, "main", "vs_5_0", flags, 0, &pBlob, &pErrorBlob);
	if (FAILED(hr))
	{
		ShowError("vertex compilation error", pErrorBlob);
		return ERROR_BAD_COMMAND;	
	}
	hr = pd3dDevice->CreateVertexShader((DWORD*)pBlob->GetBufferPointer(), pBlob->GetBufferSize(), NULL, &pVS);
	_ASSERT(SUCCEEDED(hr));
	SafeRelease(pBlob);
	SafeRelease(pErrorBlob);

//
// Pixel shader
//
	ID3D11PixelShader* pPS = NULL;
	SCOPE_EXIT(SafeRelease(pPS));
	D3D_SHADER_MACRO PSShaderMacros[] = { { "HLSL", "5.0" },{ "HLSLTOY", "1.0" }, { 0, 0 } };
	hr = D3DCompileFromFile(szArglist[1], PSShaderMacros, D3D_COMPILE_STANDARD_FILE_INCLUDE, "main", "ps_5_0", flags, 0, &pBlob, &pErrorBlob);
	if (FAILED(hr))
	{
		ShowError("pixel compilation error", pErrorBlob);
		return ERROR_BAD_COMMAND;
	}
	hr = pd3dDevice->CreatePixelShader((DWORD*)pBlob->GetBufferPointer(), pBlob->GetBufferSize(), NULL, &pPS);
	_ASSERT(SUCCEEDED(hr));
	SafeRelease(pBlob);
	SafeRelease(pErrorBlob);

//
// Uniforms buffers
//
#define HLSLTOY
#define APP_CLOUDS

#define vec2 DirectX::XMFLOAT2
#define vec3 DirectX::XMFLOAT3A
#define vec4 DirectX::XMFLOAT4
#include "../../src/uniform_buffer.h"
	main_uniform_buffer_t PSConstBuff;
	PSConstBuff.u_res.x = WIDTH;
	PSConstBuff.u_res.y = HEIGHT;
	D3D11_BUFFER_DESC uniformBuffDesc = { sizeof(PSConstBuff), D3D11_USAGE_DYNAMIC, D3D11_BIND_CONSTANT_BUFFER, D3D11_CPU_ACCESS_WRITE, 0, 0 };
	ID3D11Buffer* pUniformBuff = NULL;
	SCOPE_EXIT(SafeRelease(pUniformBuff));
	D3D11_SUBRESOURCE_DATA pData = { &PSConstBuff, 0, 0 };
	hr = pd3dDevice->CreateBuffer(&uniformBuffDesc, &pData, &pUniformBuff);
	_ASSERT(SUCCEEDED(hr));

	clouds_uniform_buffer_t clouds_settings_buff;
	D3D11_BUFFER_DESC clouds_settings_buff_descript = { sizeof(clouds_settings_buff), D3D11_USAGE_DYNAMIC, D3D11_BIND_CONSTANT_BUFFER, D3D11_CPU_ACCESS_WRITE, 0, 0 };
	ID3D11Buffer* clouds_settings_ptr = NULL;
	SCOPE_EXIT(SafeRelease(clouds_settings_ptr));
	D3D11_SUBRESOURCE_DATA clouds_settings_sub_res = { &clouds_settings_buff, 0, 0 };
	hr = pd3dDevice->CreateBuffer(&clouds_settings_buff_descript, &clouds_settings_sub_res, &clouds_settings_ptr);
	_ASSERT(SUCCEEDED(hr));

//
// Bind stuff to stages
//
	pImmediateContext->VSSetShader(pVS, NULL, 0);
	pImmediateContext->PSSetShader(pPS, NULL, 0);
	pImmediateContext->PSSetConstantBuffers(0, 1, &pUniformBuff);
	pImmediateContext->PSSetConstantBuffers(1, 1, &clouds_settings_ptr);
	pImmediateContext->PSSetShaderResources(0, 1, &pTexV);
	pImmediateContext->PSSetShaderResources(1, 1, &pNoiseTexV);
	pImmediateContext->PSSetShaderResources(2, 1, &pNoiseTexV_2);
	pImmediateContext->PSSetSamplers(0, 1, &pSampler);

//
// imgui
//
	ImGui_ImplDX11_Init(hWnd, pd3dDevice, pImmediateContext);
	SCOPE_EXIT(ImGui_ImplDX11_Shutdown());

//
// Message pump and rendering
//
	auto timerStart = std::chrono::high_resolution_clock::now();
	MSG msg = {};
	do
	{		
		if (PeekMessage(&msg, NULL, 0, 0, PM_REMOVE))
		{
			if (msg.message == WM_QUIT)
				break;
			TranslateMessage(&msg);
			DispatchMessage(&msg);
		}

		auto timerNow = std::chrono::high_resolution_clock::now();
		auto timeElapsted = std::chrono::duration_cast<std::chrono::milliseconds>(timerNow - timerStart).count();

		ImGui_ImplDX11_NewFrame();
		{
			ImGui::Text("%.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
#ifdef APP_CLOUDS
			ImGui::Spacing();
			ImGui::InputFloat3("Wind Direction", reinterpret_cast<float *>(&clouds_settings_buff.wind_dir));

			ImGui::Spacing();
			ImGui::InputFloat3("Sun Direction", reinterpret_cast<float *>(&clouds_settings_buff.sun_dir));
			ImGui::InputFloat3("Sun Colour", reinterpret_cast<float *>(&clouds_settings_buff.sun_color));
			ImGui::SliderFloat("Sun Power", &clouds_settings_buff.sun_power, 0, 12);

			ImGui::Spacing();
			ImGui::SliderInt("Ray March Steps", &clouds_settings_buff.cld_march_steps, 10, 150);
			ImGui::SliderInt("Light March Steps", &clouds_settings_buff.illum_march_steps, 1, 15);

			ImGui::Spacing();
			ImGui::InputFloat("Scaterring Coeff", &clouds_settings_buff.sigma_scattering, .05, .1);
			ImGui::SliderFloat("Coverage", &clouds_settings_buff.cld_coverage, 0., 1.);
			ImGui::SliderFloat("Thickness", &clouds_settings_buff.cld_thick, 10., 200.);
#endif
		}

		pImmediateContext->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST);
		pImmediateContext->Draw(3, 0);

		ImGui::Render();
		pSwapChain->Present(0, 0);

	// update the uniforms - use DISCARD to avoid CPU/GPU fighting for the resource
	// but this means we need to re-send everything
		PSConstBuff.u_mouse.x = sMouseButtons.x ? float(sMousePos.x) : 0.f;
		PSConstBuff.u_mouse.y = sMouseButtons.x ? float(sMousePos.y) : 0.f;
		PSConstBuff.u_time = timeElapsted / 1000.f;
		D3D11_MAPPED_SUBRESOURCE mappedResource;
		hr = pImmediateContext->Map(pUniformBuff, 0, D3D11_MAP_WRITE_DISCARD, 0, &mappedResource);
		_ASSERT(SUCCEEDED(hr));
		// must be careful not to read from data, only write
		memcpy(mappedResource.pData, &PSConstBuff, sizeof(PSConstBuff));
		pImmediateContext->Unmap(pUniformBuff, 0);

#ifdef APP_CLOUDS
		hr = pImmediateContext->Map(clouds_settings_ptr, 0, D3D11_MAP_WRITE_DISCARD, 0, &mappedResource);
		_ASSERT(SUCCEEDED(hr));
		// must be careful not to read from data, only write
		memcpy(mappedResource.pData, &clouds_settings_buff, sizeof(clouds_settings_buff));
		pImmediateContext->Unmap(clouds_settings_ptr, 0);
#endif
	} while (true);

	return ERROR_SUCCESS;
}
