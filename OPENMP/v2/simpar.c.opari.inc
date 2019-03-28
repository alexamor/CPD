#include "pomp_lib.h"

struct ompregdescr omp_rd_1 = {
  "parallelfor", "", 0, "simpar.c", 110, 110, 120, 120
};

#define POMP_DLIST_00001 shared(omp_rd_1)

struct ompregdescr omp_rd_2 = {
  "parallelfor", "", 0, "simpar.c", 130, 130, 135, 135
};

#define POMP_DLIST_00002 shared(omp_rd_2)

struct ompregdescr omp_rd_3 = {
  "parallelfor", "", 0, "simpar.c", 137, 137, 142, 142
};

#define POMP_DLIST_00003 shared(omp_rd_3)

struct ompregdescr omp_rd_4 = {
  "parallelfor", "", 0, "simpar.c", 148, 148, 154, 154
};

#define POMP_DLIST_00004 shared(omp_rd_4)

struct ompregdescr omp_rd_5 = {
  "parallelfor", "", 0, "simpar.c", 157, 157, 176, 176
};

#define POMP_DLIST_00005 shared(omp_rd_5,omp_rd_6,omp_rd_7,omp_rd_8)

struct ompregdescr omp_rd_6 = {
  "atomic", "", 0, "simpar.c", 167, 167, 168, 168
};

struct ompregdescr omp_rd_7 = {
  "atomic", "", 0, "simpar.c", 169, 169, 170, 170
};

struct ompregdescr omp_rd_8 = {
  "atomic", "", 0, "simpar.c", 173, 173, 174, 174
};

struct ompregdescr omp_rd_9 = {
  "parallelfor", "", 0, "simpar.c", 178, 178, 183, 183
};

#define POMP_DLIST_00009 shared(omp_rd_9)

struct ompregdescr omp_rd_10 = {
  "parallelfor", "", 0, "simpar.c", 186, 186, 281, 281
};

#define POMP_DLIST_00010 shared(omp_rd_10)

struct ompregdescr omp_rd_11 = {
  "parallelfor", "", 0, "simpar.c", 284, 284, 318, 318
};

#define POMP_DLIST_00011 shared(omp_rd_11)

struct ompregdescr omp_rd_12 = {
  "parallelfor", "", 0, "simpar.c", 326, 326, 338, 338
};

#define POMP_DLIST_00012 shared(omp_rd_12)

struct ompregdescr omp_rd_13 = {
  "parallelfor", "", 0, "simpar.c", 361, 361, 366, 366
};

#define POMP_DLIST_00013 shared(omp_rd_13)

struct ompregdescr omp_rd_14 = {
  "parallelfor", "", 0, "simpar.c", 370, 370, 374, 374
};

#define POMP_DLIST_00014 shared(omp_rd_14)

