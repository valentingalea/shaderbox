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
APP_EGG        | signed distance field raymarcher animation         | https://www.shadertoy.com/view/MlsGDf
APP_RAYTRACER  | simple PBR raytracer with reflection               | https://www.shadertoy.com/view/Xl2XW1
APP_SDF_AO     | test for ambient occlusion with distance fields    | https://www.shadertoy.com/view/XtBGDW
APP_2D         | old-skool demoscene tunnel effect                  |
APP_ATMOSPHERE | study of Rayleigh/Mie air scattering               | https://www.shadertoy.com/view/XtBXDz
APP_CLOUDS     | study of volumetric clouds                         | https://www.shadertoy.com/view/XtBXDw

## Requirements
* a C++11 capable compiler. Tested on
Visual Studio 2013, GCC 4.8
* any GLSL/WebGL environment. Tested on
shadertoy.com and glslsandbox.com
* any HLSL environment. Curently experimental and
only tested with fxc.exe

## Features
Over time various things were added and some
can be extracted and used separately.

* minimal raytracer framework
* minimal Cook-Torrence lighting and material setup
* signed distance field functions and operators
* 2D two bone IK solver
* a library of different noise functions
* Rayleigh/Mie atmospheric scattering solver
* a bare-bones #include expander utility to compose files
* primitive animated GIF output
