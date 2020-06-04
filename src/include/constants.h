#ifndef CONSTANTS_H
#define CONSTANTS_H

#define RESOLUTION_DEFAULT 30 /* itter */
#define RANDOM_SEED 101010

#define PI 3.1415926535897932

/* large (out-of-cache) problem sizes */

#define FFT_SIZE 1048576 /* (2^20) must be a power of two */
#define SOR_SIZE 1000    /*  NxN grid  */
#define SPARSE_SIZE_M 100000
#define SPARSE_SIZE_nz 1000000
#define LU_SIZE 1000
#define MONTECARLO_ITTER 1000000

#endif
