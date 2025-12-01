
// Cellular noise ("Worley noise") in 2D in GLSL.
// Copyright (c) Stefan Gustavson 2011-04-19. All rights reserved.
// This code is released under the conditions of the MIT license.
// See LICENSE file for details.
// https://github.com/stegu/webgl-noise

#include "common.glsl"

// Cellular noise, returning F1 and F2 in a vec2.
// Speeded up by using 2x2 search window instead of 3x3,
// at the expense of some strong pattern artifacts.
// F2 is often wrong and has sharp discontinuities.
// If you need a smooth F2, use the slower 3x3 version.
// F1 is sometimes wrong, too, but OK for most purposes.
static vec2 cellular2x2(vec2 P) {
	const float K = float(0.142857142857); // 1/7
	const float K2 = float(0.0714285714285); // K/2
	const float jitter = float(0.8); // jitter float(1.0) makes F1 wrong more often
	vec2 Pi = mod289(floor(P));
 	vec2 Pf = fract(P);
	vec4 Pfx = Pf.x + vec4(-float(0.5), -float(1.5), -float(0.5), -float(1.5));
	vec4 Pfy = Pf.y + vec4(-float(0.5), -float(0.5), -float(1.5), -float(1.5));
	vec4 p = permute(Pi.x + vec4(float(0.0), float(1.0), float(0.0), float(1.0)));
	p = permute(p + Pi.y + vec4(float(0.0), float(0.0), float(1.0), float(1.0)));
	vec4 ox = mod7(p)*K+K2;
	vec4 oy = mod7(floor(p*K))*K+K2;
	vec4 dx = Pfx + jitter*ox;
	vec4 dy = Pfy + jitter*oy;
	vec4 d = dx * dx + dy * dy; // d11, d12, d21 and d22, squared
	// Sort out the two smallest distances
#if 0
	// Cheat and pick only F1
	vec2(d.x, d.y) = min(vec2(d.x, d.y), vec2(d.z, d.w));
	d.x = min(d.x, d.y);
	return vec2(sqrt(d.x)); // F1 duplicated, F2 not computed
#else
	// Do it right and find both F1 and F2
	vec2 tmp = (d.x < d.y) ? vec2(d.x, d.y) : vec2(d.y, d.x); // Swap if smaller
	d.x = tmp.x; d.y = tmp.y;
	tmp = (d.x < d.z) ? vec2(d.x, d.z) : vec2(d.z, d.x);
	d.x = tmp.x; d.z = tmp.y;
	tmp = (d.x < d.w) ? vec2(d.x, d.w) : vec2(d.w, d.x);
	d.x = tmp.x; d.w = tmp.y;
	d.y = min(d.y, d.z);
	d.y = min(d.y, d.w);
	return sqrt(vec2(d.x, d.y));
#endif
}
