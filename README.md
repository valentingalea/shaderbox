# shader box
Various experiments with raytracing and raymarching.

Works both on C++ and GLSL using the wonderful CxxSwizzle library (https://github.com/gwiazdorrr/CxxSwizzle). Needs a C++11 compiler like Visual Studio 2013 or GCC 4.8.

The various subprojects are enabled via global defines:
* APP_EGG is an distance field ray marcher animation. (Live: https://www.shadertoy.com/view/MlsGDf)
* APP_RAYTRACER is a PBR raytracer with simple reflection and refraction.
* APP_SDF_AO is a test for ambient occlusion with signed distance fields. (Live: https://www.shadertoy.com/view/XtBGDW)
* APP_2D is an old-skool demoscene tunnel effect.
* APP_TERRAIN is an ongoing test for modeling terrain and atmosphere.
