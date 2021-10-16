# shader box
![snapshot of apps](http://i.imgur.com/JEuXZii.png)

Various experiments with raytracing and
raymarching done only in pixel (fragment) shaders.

Written in **GLSL** but supports **C++** and **HLSL** using preprocessor macro tricks.

C++ is possible due to my other library: [VML](https://github.com/valentingalea/vml)
(NOTE: needs to be placed outside this repo, in the same parent folder)

## Projects
The various submodules are enabled via global defines:

Project define | Description                                        | Live on shadertoy.com
---------------|----------------------------------------------------|-------------------------
APP_PLANET     | stylized procedural terrain and clouds             | https://www.shadertoy.com/view/ldyXRw
APP_CLOUDS     | study of volumetric clouds                         | https://www.shadertoy.com/view/XtBXDw
APP_VINYL      | vinyl turntable animation                          | https://www.shadertoy.com/view/XtG3DD
APP_EGG        | signed distance field raymarcher animation         | https://www.shadertoy.com/view/MlsGDf
APP_RAYTRACER  | simple PBR raytracer with reflection               | https://www.shadertoy.com/view/Xl2XW1
APP_ATMOSPHERE | study of Rayleigh/Mie air scattering               | https://www.shadertoy.com/view/XtBXDz

## Compiler Support
* Visual Studio C++ 2017 (15.5) +
* GCC 6.x +
* clang 3.6 +

## Environment Support
* Android - [C4droid](https://play.google.com/store/apps/details?id=com.n0n3m4.droidc&hl=en_GB)
* Web - [Shadertoy](https://www.shadertoy.com/) and [GLSL Sandbox](http://glslsandbox.com/)
* Desktop - any DirectX 11 HLSL environment; tested with included _hlsltoy_ utility 

## Features
Over time various things were added and some
can be extracted and used separately.

* minimal raytracer framework
* minimal Cook-Torrence lighting and material setup
* signed distance field functions and operators
* 2D two bone IK solver
* a library of different noise functions
* Rayleigh/Mie atmospheric scattering solver

## Utilities
The `util/` folder contains separate, independent projects:

Project   | Description
----------|-----------------------------------------------------------------------------------
inclxpnd  | a bare-bones #include expander utility to compose files
hlsltoy   | minimal one-file DX 11 framework that runs a fullscreen pixel shader
ddsvolgen | 3D noise texture generator
