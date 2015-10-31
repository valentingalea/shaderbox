#define WIDTH 1280
#define HEIGHT 720                   

#include <d3d11.h>
#include <D3Dcompiler.h>
#include <Shellapi.h>
#include <stdio.h>

static int SwapChainDesc[] = {
	WIDTH, HEIGHT, 0, 0, DXGI_FORMAT_R8G8B8A8_UNORM, 0, 0,  // w,h,refreshrate,format, scanline, scaling
	1, 0, // sampledesc (count, quality)
	DXGI_USAGE_UNORDERED_ACCESS, // usage
	1, // buffercount
	0, // hwnd
	0, 0, 0 }; // fullscreen, swap_discard, flags

static D3D11_BUFFER_DESC ConstBufferDesc = { 16, D3D11_USAGE_DEFAULT, D3D11_BIND_CONSTANT_BUFFER, 0, 0, 0 }; // 16,0,4,0,0,0

LRESULT CALLBACK WndProc(HWND hWnd, UINT uMsg, WPARAM wParam, LPARAM lParam)
{
	if (uMsg == WM_CLOSE)
		PostQuitMessage(0);
	return DefWindowProc(hWnd, uMsg, wParam, lParam);
}

int __stdcall WinMain(HINSTANCE hThisInstance, HINSTANCE hPrevInstance, LPSTR lpszArgument, int nCmdShow)
{
	LPWSTR *szArglist;
	int nArgs;
	szArglist = CommandLineToArgvW(GetCommandLineW(), &nArgs);
	if (!szArglist || nArgs < 2)
	{
		MessageBox(0, "Minimal DX11 Framework.\n\Usage:\nhlsltoy <shader file>", "Error", MB_ICONASTERISK);
		return ERROR_INVALID_COMMAND_LINE;
	}

	HRESULT hr;
	HINSTANCE hinst = GetModuleHandle(NULL);
	LPCTSTR lpszClassName = "tinyDX11";
	WNDCLASS wc = { CS_HREDRAW | CS_VREDRAW | CS_OWNDC, (WNDPROC)WndProc, 0, 0, hinst, LoadIcon(NULL, IDI_WINLOGO), LoadCursor(NULL, IDC_ARROW), 0, 0, lpszClassName };
	RegisterClass(&wc);

	HWND hWnd = CreateWindowExA(0, lpszClassName, "hlsltoy", WS_POPUPWINDOW | WS_CAPTION | WS_MINIMIZEBOX | WS_VISIBLE, 0, 0, WIDTH, HEIGHT, 0, 0, hinst, 0);
	SwapChainDesc[12] = 1; // windowed
	SwapChainDesc[11] = (int)hWnd;

	// setting up device
	ID3D11Device*               pd3dDevice;
	ID3D11DeviceContext*        pImmediateContext;
	IDXGISwapChain*             pSwapChain;
#ifdef _DEBUG	
	hr = D3D11CreateDeviceAndSwapChain(NULL, D3D_DRIVER_TYPE_HARDWARE, NULL, D3D11_CREATE_DEVICE_DEBUG, 0, 0, D3D11_SDK_VERSION, (DXGI_SWAP_CHAIN_DESC*)&SwapChainDesc[0], &pSwapChain, &pd3dDevice, NULL, &pImmediateContext);
	if (FAILED(hr)) return hr;
#else
	D3D11CreateDeviceAndSwapChain( NULL, D3D_DRIVER_TYPE_HARDWARE, NULL, 0, 0,0,D3D11_SDK_VERSION, (DXGI_SWAP_CHAIN_DESC*)&SwapChainDesc[0], &pSwapChain, &pd3dDevice, NULL, &pImmediateContext );
#endif

	// unordered access view on back buffer
	pSwapChain->GetBuffer(0, __uuidof( ID3D11Texture2D ), (LPVOID*)&SwapChainDesc[0]);	// reuse swapchain struct for temp storage of texture pointer
	pd3dDevice->CreateUnorderedAccessView((ID3D11Texture2D*)SwapChainDesc[0], NULL, (ID3D11UnorderedAccessView**)&SwapChainDesc[0]);
	pImmediateContext->CSSetUnorderedAccessViews(0, 1, (ID3D11UnorderedAccessView**)&SwapChainDesc[0], 0);

	// constant buffer
	ID3D11Buffer *pConstants;
	pd3dDevice->CreateBuffer(&ConstBufferDesc, NULL, &pConstants);
	pImmediateContext->CSSetConstantBuffers(0, 1, &pConstants);

	// compute shader
	ID3D11ComputeShader *pCS;
	ID3D10Blob* pErrorBlob = NULL;
	ID3D10Blob* pBlob = NULL;
	hr = D3DCompileFromFile(szArglist[1], NULL, D3D_COMPILE_STANDARD_FILE_INCLUDE, "cs_5_0", "cs_5_0", 0, 0, &pBlob, &pErrorBlob);
	if (FAILED(hr))
	{
		constexpr int MAX = 4096;
		char message[MAX];
		sprintf_s(message, MAX, "unknown compiler error");
		if (pErrorBlob != NULL)
			sprintf_s(message, MAX, "%s", (LPCTSTR)pErrorBlob->GetBufferPointer());
		MessageBox(hWnd, message, "Error", MB_OK);
		return 10;
	}
	else
		hr = pd3dDevice->CreateComputeShader((DWORD*)pBlob->GetBufferPointer(), pBlob->GetBufferSize(), NULL, &pCS);
	void *pp = pBlob->GetBufferPointer(); int ps = pBlob->GetBufferSize();

	if (pErrorBlob != 0) pErrorBlob->Release();
	if (pBlob != 0) pBlob->Release();
	pImmediateContext->CSSetShader(pCS, NULL, 0);

	MSG msg;
	do
	{		
		if (PeekMessage(&msg, NULL, 0, 0, PM_REMOVE))
		{
			if (msg.message == WM_QUIT)
				break;
			TranslateMessage(&msg);
			DispatchMessage(&msg);
		}

		SwapChainDesc[0] = GetTickCount();
		pImmediateContext->UpdateSubresource(pConstants, 0, 0, &SwapChainDesc[0], 4, 4);
		pImmediateContext->Dispatch(WIDTH / 16, HEIGHT / 16, 1);
		pSwapChain->Present(0, 0);
	} while (true);

	LocalFree(szArglist);
	ExitProcess(0);
	return ERROR_SUCCESS;
}
