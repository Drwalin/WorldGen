#pragma once

#define SAFE_COND_GRID(COND, EXPR) \
	if constexpr (safe) { \
		if(COND) { \
			EXPR; \
		} \
	} else { \
		EXPR; \
	}

#define NEIGHBOURS(VAR_NAME, X, Y) \
	int VAR_NAME[4] = { \
		Neighbour<safe, 0>(X, Y), \
		Neighbour<safe, 1>(X, Y), \
		Neighbour<safe, 2>(X, Y), \
		Neighbour<safe, 3>(X, Y) }

#define FOR_EACH_DIR(CODE) \
{ \
	{ \
		constexpr int DIR = 0; \
		constexpr int R_DIR = 2; \
		(void)DIR; (void)R_DIR; \
		CODE; \
	} \
	{ \
		constexpr int DIR = 1; \
		constexpr int R_DIR = 3; \
		(void)DIR; (void)R_DIR; \
		CODE; \
	} \
	{ \
		constexpr int DIR = 2; \
		constexpr int R_DIR = 0; \
		(void)DIR; (void)R_DIR; \
		CODE; \
	} \
	{ \
		constexpr int DIR = 3; \
		constexpr int R_DIR = 1; \
		(void)DIR; (void)R_DIR; \
		CODE; \
	} \
}

#define FOR_EACH_DIR_SAFE_COND(COND, CODE) \
{ \
	{ \
		constexpr int DIR = 0; \
		constexpr int R_DIR = 2; \
		(void)DIR; (void)R_DIR; \
		SAFE_COND_GRID(COND, CODE); \
	} \
	{ \
		constexpr int DIR = 1; \
		constexpr int R_DIR = 3; \
		(void)DIR; (void)R_DIR; \
		SAFE_COND_GRID(COND, CODE); \
	} \
	{ \
		constexpr int DIR = 2; \
		constexpr int R_DIR = 0; \
		(void)DIR; (void)R_DIR; \
		SAFE_COND_GRID(COND, CODE); \
	} \
	{ \
		constexpr int DIR = 3; \
		constexpr int R_DIR = 1; \
		(void)DIR; (void)R_DIR; \
		SAFE_COND_GRID(COND, CODE); \
	} \
}

#define REV_DIR(DIR) \
	((DIR + 2) & 3)

