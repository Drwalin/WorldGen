#include <random>

#include "../include/worldgen/Hash.hpp"

namespace wg
{
namespace hx
{
static uint8_t lookUpTable[256][256];
static uint8_t fill = []() -> uint8_t {
	std::mt19937_64 mt(12345);
	for (int i = 0; i < 256; ++i) {
		for (int j=0; j<256; ++j) {
			lookUpTable[i][j] = mt();
		}
	}
	return 0;
}();

static inline uint64_t ror64(uint64_t v, int r)
{
	return (v >> r) | (v << (64 - r));
}
static inline float ToFloat(uint64_t h)
{
	h = mxrmx(ror64(h, 7));
	constexpr uint64_t prime = 24432403ull;
	uint64_t r = h % prime;
	return ((double)r) / (double)(prime - 1);
}

float Hash(glm::ivec3 p, int seed)
{
	uint64_t h1 = mxrmx(mxrmx((uint64_t)p.x) | mxrmx(ror64((uint64_t)p.y, 32)));
	uint64_t h2 = mxrmx(mxrmx((uint64_t)p.y) | mxrmx(ror64((uint64_t)p.z, 32)));
	uint64_t h3 = mxrmx(mxrmx((uint64_t)p.z) | mxrmx(ror64((uint64_t)p.x, 32)));

	uint64_t h = h1 ^ ror64(h2, 17) ^ ror64(h3, 29) ^ mxrmx(seed);
	return ToFloat(h);
}

float Hash(glm::ivec2 p, int seed)
{
	uint64_t v = 0;
	v = lookUpTable[(seed>>16) & 255][(seed >> 24) & 255];
	v = lookUpTable[(seed) & 255][(v + (seed >> 8)) & 255];
	for (int i = 3; i >= 0; --i) {
		v = lookUpTable[(v ^ (p.x>>(i<<3))) & 255][(v + (p.y>>(i<<3))) & 255];
	}
	return v / 255.0f;
	
	
	uint64_t h = mxrmx((mxrmx((uint64_t)p.x) | mxrmx(ror64((uint64_t)p.y, 32))) ^ seed);
	return ToFloat(h);
}

glm::vec2 Hash2(glm::ivec2 p, int seed)
{
	uint64_t h = mxrmx(((uint64_t)p.x) | (((uint64_t)p.y) << 32)) ^ mxrmx(seed);
	
	return {ToFloat(h), ToFloat(ror64(h, 32))};
}

// x c2 mul 56 32 xrr c3 mul 23 xsr
uint64_t mxrmx(uint64_t x)
{
	constexpr int HASH = 0;
	if (HASH != 5) {
		x ^= 0xcbf29ce484222325ull;
		x *= 0x00000100000001b3ull;
		x = ror64(x, 27);
		x *= 0x00000100000001b3ull;
	}
	switch (HASH) {
	case 0:
		x *= 0xbf58476d1ce4e5b9ull;
		x ^= x >> 32;
		x *= 0x94d049bb133111ebull;
		x ^= x >> 32;
		x *= 0x94d049bb133111ebull;
		return x;

	case 1:
		x *= 0x94d049bb133111ebull;
		x ^= ror64(x, 56) ^ ror64(x, 32);
		x *= 0xff51afd7ed558ccdull;
		x ^= x >> 23;
		return x;
	case 2:

		x *= 0xbf58476d1ce4e5b9ull;
		x ^= x >> 56;
		x *= 0x94d049bb133111ebull;
		return x;

	case 3:
		x ^= x >> 33;
		x *= 0xff51afd7ed558ccdull;
		x ^= x >> 33;
		x *= 0xc4ceb9fe1a85ec53ull;
		x ^= x >> 33;
		return x;

	case 4: {
		uint64_t hash = 0;
		uint64_t b0 = (x & 255);
		uint64_t b1 = (x & 65280) >> 8;
		uint64_t b2 = (x & 16711680) >> 16;
		uint64_t b3 = (x & -16777216) >> 24;
		hash += b0;
		hash += (hash << 10);
		hash ^= (hash >> 6);
		hash += b1;
		hash += (hash << 10);
		hash ^= (hash >> 6);
		hash += b2;
		hash += (hash << 10);
		hash ^= (hash >> 6);
		hash += b3;
		hash += (hash << 10);
		hash ^= (hash >> 6);
		hash += (hash << 3);
		hash ^= (hash >> 11);
		hash += (hash << 15);
		return hash;
	}

	case 5: {
		uint64_t hash = 0;
		uint8_t *b = (uint8_t *)&x;
		uint8_t *d = (uint8_t *)&hash;

		for (int i = 0; i < 8; ++i) {
			uint8_t v = i;
			for (int j = 1; j < 8; ++j) {
				v = lookUpTable[v][b[(j + i)&7]];
			}
			d[i] = v;
		}
		return hash;
	}
	}
}
} // namespace hx
} // namespace wg
