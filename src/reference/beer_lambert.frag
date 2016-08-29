// https://www.shadertoy.com/view/Mlc3Ds

float BeerLambert(float sigma_t, float d)
{
    return exp(-sigma_t * d);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    // compute length of medium from mouse X
    float d = iMouse.x / iResolution.x * 5.0;
    // compute extinction factor from mouse Y
	float sigma_t = iMouse.y / iResolution.y * 20.0;
    
    // default input (for shadertoy preview)
    if (iMouse.xy == vec2(0.0)) {
        d = 2.0;
        sigma_t = 5.0;
    }
    
    // how far the current pixel has traveled in the medium
    float x = fragCoord.x / iResolution.x * d;
    // compute the extinction function for this pixel's distance in the medium
    float T = BeerLambert(sigma_t, x);
    
    // convert to a color
    vec4 color = vec4(vec3(T),1.0);
    // gamma-correct it
    color.rgb = pow(color.rgb, vec3(1.0/2.2));
    
    // display color
    fragColor = color;
}