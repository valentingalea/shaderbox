#include "def.h"

// Created by Reinder Nijhoff 2015
// @reindernijhoff

void mainImage(_out(vec4) f, _in(vec2) w)
{
	vec3 d = (vec3(w.xy, 1.) / iResolution.x)  - .5;
	vec3 p, c;
	vec3 g = d, o = d;
	o.z += fract(sin(iGlobalTime)) * 4.;

	for (float i = .0; i < 9.; i += .01)
	{
		p = (c = o += d * i * .05) * .3;
		if (cos(p.z) - abs(sin(p.x * .7 + cos(p.z))) > ++p.y)
		{
			g = mix((3. + p.y) * vec3(.6, .3, 0), d, i / 9.);
			break;
		}
	}
	f.xyz = g;
}