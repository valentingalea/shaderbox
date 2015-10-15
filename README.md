# shader box
Various experiments with raytracing and raymarching.

Works both on C++ and GLSL using the wonderful CxxSwizzle library (https://github.com/gwiazdorrr/CxxSwizzle). Needs a C++11 compiler like Visual Studio 2013 or GCC 4.8.

The various subprojects are enabled via global defines:

Project define | Description                                            | Live on shadertoy.com
---------------|--------------------------------------------------------|-------------------
APP_EGG        | an distance field ray marcher animation                | https://www.shadertoy.com/view/MlsGDf
APP_RAYTRACER  | PBR raytracer with simple reflection and refraction    | https://www.shadertoy.com/view/Xl2XW1
APP_SDF_AO     | test for ambient occlusion with signed distance fields | https://www.shadertoy.com/view/XtBGDW
APP_2D         | old-skool demoscene tunnel effect                      |
APP_ATMOSPHERE | study of Rayleigh/Mie air scattering                   | https://www.shadertoy.com/view/XtBXDz
APP_CLOUDS     | study of volumetric clouds                             | https://www.shadertoy.com/view/XtBXDw

