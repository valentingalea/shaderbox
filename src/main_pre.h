// CxxSwizzle
// Copyright (c) 2013, Piotr Gwiazdowski <gwiazdorrr+github at gmail.com>

#include <swizzle/glsl/naive/vector.h>
#include <swizzle/glsl/naive/matrix.h>
#include <swizzle/glsl/texture_functions.h>

typedef swizzle::glsl::naive::vector< int, 2 > ivec2;
typedef swizzle::glsl::naive::vector< float, 2 > vec2;
typedef swizzle::glsl::naive::vector< float, 3 > vec3;
typedef swizzle::glsl::naive::vector< float, 4 > vec4;

static_assert(sizeof(vec2) == sizeof(float[2]), "Too big");
static_assert(sizeof(vec3) == sizeof(float[3]), "Too big");
static_assert(sizeof(vec4) == sizeof(float[4]), "Too big");

typedef swizzle::glsl::naive::matrix< swizzle::glsl::naive::vector, float, 2, 2> mat2;
typedef swizzle::glsl::naive::matrix< swizzle::glsl::naive::vector, float, 3, 3> mat3;
typedef swizzle::glsl::naive::matrix< swizzle::glsl::naive::vector, float, 4, 4> mat4;

//! A really, really simplistic sampler using SDLImage
struct SDL_Surface;
class sampler2D : public swizzle::glsl::texture_functions::tag
{
public:
    enum WrapMode
    {
        Clamp,
        Repeat,
        MirrorRepeat
    };

    typedef vec2 tex_coord_type;

    sampler2D(const char* path, WrapMode wrapMode);
    ~sampler2D();
    vec4 sample(vec2 coord);

private:
    SDL_Surface *m_image;
    WrapMode m_wrapMode;

    // do not allow copies to be made
    sampler2D(const sampler2D&);
    sampler2D& operator=(const sampler2D&);
};

// this where the magic happens...
namespace glsl_sandbox
{
    // a nested namespace used when redefining 'inout' and 'out' keywords
    namespace ref
    {
        typedef ::vec2& vec2;
        typedef ::vec3& vec3;
        typedef ::vec4& vec4;
    }

    #include <swizzle/glsl/vector_functions.h>

    // constants shaders are using
    float time;
    vec2 mouse;
    vec2 resolution;

    // constants some shaders from shader toy are using
    vec2& iResolution = resolution;
    float& iGlobalTime = time;
    vec2& iMouse = mouse;

    sampler2D diffuse("diffuse.png", sampler2D::Repeat);
    sampler2D specular("specular.png", sampler2D::Repeat);

    struct fragment_shader
    {
        vec2 gl_FragCoord;
        vec4 gl_FragColor;
        void operator()(void);
    };

    // change meaning of glsl keywords to match sandbox
    #define uniform extern
    #define in
    #define out ref::
    #define inout ref::
    #define main fragment_shader::operator()

    #pragma warning(push)
    #pragma warning(disable: 4244) // disable return implicit conversion warning
    #pragma warning(disable: 4305) // disable truncation warning
    
