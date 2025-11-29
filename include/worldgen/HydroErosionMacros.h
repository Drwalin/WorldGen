
#define COND_GRID(COND, EXPR)                                                  \
	if (COND) {                                                                \
		EXPR                                                                   \
	}

#define NEIGHBOURS(VAR_NAME, X, Y)                                             \
	int VAR_NAME[4] = {Neighbour(X, Y, 0), Neighbour(X, Y, 1),                 \
					   Neighbour(X, Y, 2), Neighbour(X, Y, 3)}

#define NEIGHBOURS_CORNERS(VAR_NAME, X, Y)                                     \
	int VAR_NAME[4] = {                                                        \
		At(X - 1, Y - 1),                                                      \
		At(X + 1, Y - 1),                                                      \
		At(X - 1, Y + 1),                                                      \
		At(X + 1, Y + 1),                                                      \
	}

#define FOR_EACH_DIR(CODE)                                                     \
	{                                                                          \
		for (int DIR = 0; DIR < 4; ++DIR) {                                    \
			int R_DIR = (DIR + 2) & 3;                                         \
			IGNORE_UNUSED(R_DIR)                                               \
			CODE                                                               \
		}                                                                      \
	}

#define FOR_EACH_DIR_COND(COND, CODE) FOR_EACH_DIR(COND_GRID(COND, CODE))

#define REV_DIR(DIR) ((DIR + 2) & 3)

#define SWAP(A, B)                                                             \
	{                                                                          \
		float tmp = A;                                                         \
		A = B;                                                                 \
		B = tmp;                                                               \
	}
