# shader box
![snapshot of apps](http://valentingalea.github.io/shaderbox/apps.png)

Various experiments with raytracing and
raymarching done only in pixel (fragment) shaders.

Written in GLSL but supports C++ and HLSL
using preprocessor macro tricks.
The C++11 work is made possible by the
wonderful CxxSwizzle library.
(https://github.com/gwiazdorrr/CxxSwizzle)

## Projects
The various submodules are enabled via global defines:

Project define | Description                                        | Live on shadertoy.com
---------------|----------------------------------------------------|-------------------------
APP_PLANET     | stylized procedural terrain and clouds             | https://www.shadertoy.com/view/ldyXRw
APP_ATMOSPHERE | study of Rayleigh/Mie air scattering               | https://www.shadertoy.com/view/XtBXDz
APP_CLOUDS     | study of volumetric clouds                         | https://www.shadertoy.com/view/XtBXDw
APP_EGG        | signed distance field raymarcher animation         | https://www.shadertoy.com/view/MlsGDf
APP_RAYTRACER  | simple PBR raytracer with reflection               | https://www.shadertoy.com/view/Xl2XW1
APP_SDF_AO     | test for ambient occlusion with distance fields    | https://www.shadertoy.com/view/XtBGDW

## Support
* Visual Studio C++ 2015. Has a dependency on SDL 1.2 (included).
* [C4droid](https://play.google.com/store/apps/details?id=com.n0n3m4.droidc&hl=en_GB) on Android.
* any GLSL ES environment. Tested on [Shadertoy](https://www.shadertoy.com/) and [GLSL Sandbox](http://glslsandbox.com/).
* any DirectX 11 HLSL environment. Tested with included hlsltoy utility. Has a dependency on _DirectXTex_ library (external repository reference).

## Features
Over time various things were added and some
can be extracted and used separately.

* minimal raytracer framework
* minimal Cook-Torrence lighting and material setup
* signed distance field functions and operators
* 2D two bone IK solver
* a library of different noise functions
* Rayleigh/Mie atmospheric scattering solver
* primitive animated GIF output

## Utilities
The util/ folder contains separate, independent projects:

Project   | Description
----------|-----------------------------------------------------------------------------------
inclxpnd  | a bare-bones #include expander utility to compose files
hlsltoy   | minimal one-file DX 11 framework that runs a fullscreen pixel shader
ddsvolgen | 3D noise texture generator
