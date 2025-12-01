
// Cellular noise ("Worley noise") in 2D in GLSL.
// Copyright (c) Stefan Gustavson 2011-04-19. All rights reserved.
// This code is released under the conditions of the MIT license.
// See LICENSE file for details.
// https://github.com/stegu/webgl-noise

#include "common.glsl"

// Cellular noise, returning F1 and F2 in a vec2.
// Standard 3x3 search window for good F1 and F2 values
static vec2 cellular(vec2 P)
{
	const float K = float(0.142857142857);  // 1/7
	const float Ko = float(0.428571428571); // 3/7
	const float jitter = float(1.0);		  // Less gives more regular pattern
	vec2 Pi = mod289(floor(P));
	vec2 Pf = fract(P);
	vec3 oi = vec3(-float(1.0), float(0.0), float(1.0));
	vec3 of = vec3(-float(0.5), float(0.5), float(1.5));
	vec3 px = permute(Pi.x + oi);
	vec3 p = permute(px.x + Pi.y + oi); // p11, p12, p13
	vec3 ox = fract(p * K) - Ko;
	vec3 oy = mod7(floor(p * K)) * K - Ko;
	vec3 dx = Pf.x + float(0.5) + jitter * ox;
	vec3 dy = Pf.y - of + jitter * oy;
	vec3 d1 = dx * dx + dy * dy;   // d11, d12 and d13, squared
	p = permute(px.y + Pi.y + oi); // p21, p22, p23
	ox = fract(p * K) - Ko;
	oy = mod7(floor(p * K)) * K - Ko;
	dx = Pf.x - float(0.5) + jitter * ox;
	dy = Pf.y - of + jitter * oy;
	vec3 d2 = dx * dx + dy * dy;   // d21, d22 and d23, squared
	p = permute(px.z + Pi.y + oi); // p31, p32, p33
	ox = fract(p * K) - Ko;
	oy = mod7(floor(p * K)) * K - Ko;
	dx = Pf.x - float(1.5) + jitter * ox;
	dy = Pf.y - of + jitter * oy;
	vec3 d3 = dx * dx + dy * dy; // d31, d32 and d33, squared
	// Sort out the two smallest distances (F1, F2)
	vec3 d1a = min(d1, d2);
	d2 = max(d1, d2);					   // Swap to keep candidates for F2
	d2 = min(d2, d3);					   // neither F1 nor F2 are now in d3
	d1 = min(d1a, d2);					   // F1 is now in d1
	d2 = max(d1a, d2);					   // Swap to keep candidates for F2
	vec2 tmp = (d1.x < d1.y) ? vec2(d1.x, d1.y) : vec2(d1.y, d1.x); // Swap if smaller
	d1.x = tmp.x; d1.y = tmp.y;
	tmp = (d1.x < d1.z) ? vec2(d1.x, d1.z) : vec2(d1.z, d1.x); // F1 is in d1.x
	d1.x = tmp.x; d1.z = tmp.y;
	tmp = min(vec2(d1.y, d1.z), vec2(d2.y, d2.z));			   // F2 is now not in vec2(d2.y, d2.z)
	d1.y = tmp.x; d1.z = tmp.y;
	d1.y = min(d1.y, d1.z);				   // nor in  d1.z
	d1.y = min(d1.y, d2.x);				   // F2 is in d1.y, we're done.
	return sqrt(vec2(d1.x, d1.y));
}
