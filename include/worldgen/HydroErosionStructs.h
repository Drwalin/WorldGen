
#if !defined HYDRO_EROSION_STRUCTS_DEFINED
#define HYDRO_EROSION_STRUCTS_DEFINED
struct Flux {
	float f[4];
};

struct Velocity {
	float x;
	float y;
};

struct GroundLayers {
	float layers[2];
};
#endif

