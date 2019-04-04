#include "pomp_lib.h"

struct ompregdescr omp_rd_1 = {
  "parallelfor", "", 0, "simpar-omp.c", 111, 111, 121, 121
};

#define POMP_DLIST_00001 shared(omp_rd_1)

struct ompregdescr omp_rd_2 = {
  "parallelfor", "", 0, "simpar-omp.c", 128, 128, 134, 134
};

#define POMP_DLIST_00002 shared(omp_rd_2)

struct ompregdescr omp_rd_3 = {
  "parallelfor", "", 0, "simpar-omp.c", 137, 137, 156, 156
};

#define POMP_DLIST_00003 shared(omp_rd_3,omp_rd_4,omp_rd_5,omp_rd_6)

struct ompregdescr omp_rd_4 = {
  "atomic", "", 0, "simpar-omp.c", 147, 147, 148, 148
};

struct ompregdescr omp_rd_5 = {
  "atomic", "", 0, "simpar-omp.c", 149, 149, 150, 150
};

struct ompregdescr omp_rd_6 = {
  "atomic", "", 0, "simpar-omp.c", 153, 153, 154, 154
};

struct ompregdescr omp_rd_7 = {
  "parallelfor", "", 0, "simpar-omp.c", 158, 158, 163, 163
};

#define POMP_DLIST_00007 shared(omp_rd_7)

struct ompregdescr omp_rd_8 = {
  "parallelfor", "", 0, "simpar-omp.c", 166, 166, 270, 270
};

#define POMP_DLIST_00008 shared(omp_rd_8)

struct ompregdescr omp_rd_9 = {
  "parallelfor", "", 0, "simpar-omp.c", 273, 273, 307, 307
};

#define POMP_DLIST_00009 shared(omp_rd_9)

struct ompregdescr omp_rd_10 = {
  "parallelfor", "", 0, "simpar-omp.c", 315, 315, 327, 327
};

#define POMP_DLIST_00010 shared(omp_rd_10)

struct ompregdescr omp_rd_11 = {
  "parallelfor", "", 0, "simpar-omp.c", 352, 352, 355, 355
};

#define POMP_DLIST_00011 shared(omp_rd_11)

