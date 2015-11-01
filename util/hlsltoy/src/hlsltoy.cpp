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

#define SafeRelease(T) if (T) { T->Release(); }

#include <d3d11.h>
#include <D3Dcompiler.h>
#include <DirectXMath.h>
#include <CRTDBG.H>
#include <cstdio>
#include <chrono>

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
constexpr int WIDTH = 800;
constexpr int HEIGHT = 600;

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

LRESULT CALLBACK WndProc(HWND hWnd, UINT uMsg, WPARAM wParam, LPARAM lParam)
{
	if (uMsg == WM_CLOSE)
		PostQuitMessage(0);
	return DefWindowProc(hWnd, uMsg, wParam, lParam);
}

int __stdcall WinMain(HINSTANCE hThisInstance, HINSTANCE hPrevInstance, LPSTR lpszArgument, int nCmdShow)
{
	LPWSTR *szArglist = NULL;
	SCOPE_EXIT(if (szArglist) LocalFree(szArglist));
	int nArgs = 0;
	szArglist = CommandLineToArgvW(GetCommandLineW(), &nArgs);
	if (!szArglist || nArgs < 2)
	{
		ShowError("Minimal DX11 Framework.\n\Usage:\nhlsltoy <shader file>");
		return ERROR_INVALID_COMMAND_LINE;
	}

	HRESULT hr = NULL;
	HINSTANCE hinst = GetModuleHandle(NULL);
	WNDCLASS wc = { CS_HREDRAW | CS_VREDRAW | CS_OWNDC, (WNDPROC)WndProc, 0, 0, hinst, LoadIcon(NULL, IDI_WINLOGO), LoadCursor(NULL, IDC_ARROW), 0, 0, lpszClassName };
	RegisterClass(&wc);
	SCOPE_EXIT(UnregisterClass(lpszClassName, hinst));

	HWND hWnd = CreateWindowExA(0, lpszClassName, lpszAppName, WS_POPUPWINDOW | WS_CAPTION | WS_MINIMIZEBOX | WS_VISIBLE, CW_USEDEFAULT, CW_USEDEFAULT, WIDTH, HEIGHT, 0, 0, hinst, 0);
	SwapChainDesc.OutputWindow = hWnd;

	// setting up device
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

	ID3D11Texture2D* pBackBuffer = NULL;
	SCOPE_EXIT(SafeRelease(pBackBuffer));
	hr = pSwapChain->GetBuffer(0, __uuidof(ID3D11Texture2D), (LPVOID*)&pBackBuffer);
	_ASSERT(SUCCEEDED(hr));

	ID3D11RenderTargetView* pRenderTargetView = NULL;
	SCOPE_EXIT(SafeRelease(pRenderTargetView));
	hr = pd3dDevice->CreateRenderTargetView(pBackBuffer, NULL, &pRenderTargetView);
	_ASSERT(SUCCEEDED(hr));

	pImmediateContext->OMSetRenderTargets(1, &pRenderTargetView, NULL);

	D3D11_VIEWPORT vp = { 0, 0, WIDTH, HEIGHT, 0., 1. };
	pImmediateContext->RSSetViewports(1, &vp);

	D3D11_RASTERIZER_DESC rasterizerDesc = { D3D11_FILL_SOLID, D3D11_CULL_NONE, FALSE, 0, 0., 0., FALSE, FALSE, FALSE, FALSE };
	ID3D11RasterizerState* pd3dRasterizerState = NULL;
	SCOPE_EXIT(SafeRelease(pd3dRasterizerState));
	hr = pd3dDevice->CreateRasterizerState(&rasterizerDesc, &pd3dRasterizerState);
	_ASSERT(SUCCEEDED(hr));

	ID3D10Blob* pErrorBlob = NULL;
	ID3D10Blob* pBlob = NULL;
	SCOPE_EXIT(SafeRelease(pErrorBlob));
	SCOPE_EXIT(SafeRelease(pBlob));
	UINT flags = D3DCOMPILE_ENABLE_STRICTNESS | D3DCOMPILE_IEEE_STRICTNESS | D3DCOMPILE_OPTIMIZATION_LEVEL3 | D3DCOMPILE_WARNINGS_ARE_ERRORS;

	// vertex shader
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

	// pixel shader
	ID3D11PixelShader* pPS = NULL;
	SCOPE_EXIT(SafeRelease(pPS));
	//D3D_SHADER_MACRO PSShaderMacros[1] = { "__HLSL", "1" };
	hr = D3DCompileFromFile(szArglist[1], NULL, D3D_COMPILE_STANDARD_FILE_INCLUDE, "main", "ps_5_0", 0, 0, &pBlob, &pErrorBlob);
	if (FAILED(hr))
	{
		ShowError("pixel compilation error", pErrorBlob);
		return ERROR_BAD_COMMAND;
	}
	hr = pd3dDevice->CreatePixelShader((DWORD*)pBlob->GetBufferPointer(), pBlob->GetBufferSize(), NULL, &pPS);
	_ASSERT(SUCCEEDED(hr));
	SafeRelease(pBlob);
	SafeRelease(pErrorBlob);

	__declspec(align(16)) struct PS_CONSTANT_BUFFER
	{
		DirectX::XMFLOAT2 resolution;
		float time;
		DirectX::XMFLOAT2 mouse;
	};
	PS_CONSTANT_BUFFER PSConstBuff = { {float(WIDTH), float(HEIGHT)}, 0, {0, 0} };
	D3D11_BUFFER_DESC uniformBuffDesc = { sizeof(PS_CONSTANT_BUFFER), D3D11_USAGE_DYNAMIC, D3D11_BIND_CONSTANT_BUFFER, D3D11_CPU_ACCESS_WRITE, 0, 0 };
	ID3D11Buffer* pUniformBuff = NULL;
	SCOPE_EXIT(SafeRelease(pUniformBuff));
	D3D11_SUBRESOURCE_DATA pData = { &PSConstBuff, 0, 0 };
	hr = pd3dDevice->CreateBuffer(&uniformBuffDesc, &pData, &pUniformBuff);
	_ASSERT(SUCCEEDED(hr));

	pImmediateContext->VSSetShader(pVS, NULL, 0);
	pImmediateContext->PSSetShader(pPS, NULL, 0);
	pImmediateContext->PSSetConstantBuffers(0, 1, &pUniformBuff);

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

#if 0 // TODO: doesn't seem to work - CPU/GPU trample over each other
		// update the uniforms - potential bottleneck here as CPU/GPU fight for the resource
		D3D11_MAPPED_SUBRESOURCE mappedResource;
		pImmediateContext->Map(pUniformBuff, 0, D3D11_MAP_WRITE_DISCARD, 0, &mappedResource);
		// must be careful not to read from data, only write
		volatile PS_CONSTANT_BUFFER *pBuff = (PS_CONSTANT_BUFFER *)mappedResource.pData;
		pBuff->time = (float)(timeElapsted / 1000.f);
		pImmediateContext->Unmap(pUniformBuff, 0);
#endif

		pImmediateContext->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_TRIANGLELIST);
		pImmediateContext->Draw(3, 0);

		pSwapChain->Present(0, 0);
	} while (true);
	
	return ERROR_SUCCESS;
}
