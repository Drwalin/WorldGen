//
// Description : Array and textureless GLSL 2D/3D/4D simplex 
//               noise functions.
//      Author : Ian McEwan, Ashima Arts.
//  Maintainer : stegu
//     Lastmod : 20110822 (ijm)
//     License : Copyright (C) 2011 Ashima Arts. All rights reserved.
//               Distributed under the MIT License. See LICENSE file.
//               https://github.com/ashima/webgl-noise
//               https://github.com/stegu/webgl-noise
// 

#include "common.glsl"

vec4 grad4(float j, vec4 ip)
  {
  const vec4 ones = vec4(float(1.0), float(1.0), float(1.0), -float(1.0));
  vec4 p,s;

  vec3 tmp = floor( fract (vec3(j) * vec3(ip.x, ip.y, ip.z)) * float(7.0)) * ip.z - float(1.0);
  p.x = tmp.x; p.y = tmp.y; p.z = tmp.z;
  p.w = float(1.5) - dot(abs(vec3(p.x, p.y, p.z)), vec3(ones.x, ones.y, ones.z));
  s = vec4(lessThan(p, vec4(float(0.0))));
  tmp = vec3(p.x, p.y, p.z) + (vec3(s.x, s.y, s.z)*float(2.0) - float(1.0)) * vec3(s.w, s.w, s.w); 
  p.x = tmp.x; p.y = tmp.y; p.z = tmp.z;

  return p;
  }
						
// (sqrt(5) - 1)/4 = F4, used once below
#define F4 float(0.309016994374947451)

float snoise(vec4 v)
  {
  const vec4  C = vec4( float(0.138196601125011),  // (5 - sqrt(5))/20  G4
                        float(0.276393202250021),  // 2 * G4
                        float(0.414589803375032),  // 3 * G4
                       -float(0.447213595499958)); // -1 + 4 * G4

// First corner
  vec4 i  = floor(v + dot(v, vec4(F4)) );
  vec4 x0 = v -   i + dot(i, vec4(C.x, C.x, C.x, C.x));

// Other corners

// Rank sorting originally contributed by Bill Licea-Kane, AMD (formerly ATI)
  vec4 i0;
  vec3 isX = step( vec3(x0.y, x0.z, x0.w), vec3(x0.x, x0.x, x0.x) );
  vec3 isYZ = step( vec3(x0.z, x0.w, x0.w), vec3(x0.y, x0.y, x0.z) );
//  i0.x = dot( isX, vec3( float(1.0) ) );
  i0.x = isX.x + isX.y + isX.z;
  vec3 tmp = float(1.0) - isX;
  i0.y = tmp.x; i0.z = tmp.y; i0.w = tmp.z;
//  i0.y += dot( isYZ.xy, vec2( float(1.0) ) );
  i0.y += isYZ.x + isYZ.y;
  i0.z += float(1.0) - isYZ.x;
  i0.w += float(1.0) - isYZ.y;
  
  
  
  i0.z += isYZ.z;
  i0.w += float(1.0) - isYZ.z;

  // i0 now contains the unique values 0,1,2,3 in each channel
  vec4 i3 = clamp( i0, float(0.0), float(1.0) );
  vec4 i2 = clamp( i0-float(1.0), float(0.0), float(1.0) );
  vec4 i1 = clamp( i0-float(2.0), float(0.0), float(1.0) );

  //  x0 = x0 - float(0.0) + float(0.0) * vec4(C.x, C.x, C.x, C.x)
  //  x1 = x0 - i1  + float(1.0) * vec4(C.x, C.x, C.x, C.x)
  //  x2 = x0 - i2  + float(2.0) * vec4(C.x, C.x, C.x, C.x)
  //  x3 = x0 - i3  + float(3.0) * vec4(C.x, C.x, C.x, C.x)
  //  x4 = x0 - float(1.0) + float(4.0) * vec4(C.x, C.x, C.x, C.x)
  vec4 x1 = x0 - i1 + vec4(C.x, C.x, C.x, C.x);
  vec4 x2 = x0 - i2 + vec4(C.y, C.y, C.y, C.y);
  vec4 x3 = x0 - i3 + vec4(C.z, C.z, C.z, C.z);
  vec4 x4 = x0 + vec4(C.w, C.w, C.w, C.w);

// Permutations
  i = mod289(i); 
  float j0 = permute( permute( permute( permute(i.w) + i.z) + i.y) + i.x);
  vec4 j1 = permute( permute( permute( permute (
             i.w + vec4(i1.w, i2.w, i3.w, float(1.0) ))
           + i.z + vec4(i1.z, i2.z, i3.z, float(1.0) ))
           + i.y + vec4(i1.y, i2.y, i3.y, float(1.0) ))
           + i.x + vec4(i1.x, i2.x, i3.x, float(1.0) ));

// Gradients: 7x7x6 points over a cube, mapped onto a 4-cross polytope
// 7*7*6 = 294, which is close to the ring size 17*17 = 289.
  vec4 ip = vec4(float(1.0)/float(294.0), float(1.0)/float(49.0), float(1.0)/float(7.0), float(0.0)) ;

  vec4 p0 = grad4(j0,   ip);
  vec4 p1 = grad4(j1.x, ip);
  vec4 p2 = grad4(j1.y, ip);
  vec4 p3 = grad4(j1.z, ip);
  vec4 p4 = grad4(j1.w, ip);

// Normalise gradients
  vec4 norm = taylorInvSqrt(vec4(dot(p0,p0), dot(p1,p1), dot(p2, p2), dot(p3,p3)));
  p0 *= norm.x;
  p1 *= norm.y;
  p2 *= norm.z;
  p3 *= norm.w;
  p4 *= taylorInvSqrt(dot(p4,p4));

// Mix contributions from the five corners
  vec3 m0 = max(float(0.57) - vec3(dot(x0,x0), dot(x1,x1), dot(x2,x2)), float(0.0));
  vec2 m1 = max(float(0.57) - vec2(dot(x3,x3), dot(x4,x4)            ), float(0.0));
  m0 = m0 * m0;
  m1 = m1 * m1;
  return float(60.1) * ( dot(m0*m0, vec3( dot( p0, x0 ), dot( p1, x1 ), dot( p2, x2 )))
               + dot(m1*m1, vec2( dot( p3, x3 ), dot( p4, x4 ) ) ) ) ;

  }
