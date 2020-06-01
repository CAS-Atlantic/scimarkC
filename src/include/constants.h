#ifndef CONSTANTS_H_
#define CONSTANTS_H_

const double RESOLUTION_DEFAULT = 30; /* itter */
const unsigned int RANDOM_SEED = 101010;

/* large (out-of-cache) problem sizes */

const unsigned int FFT_SIZE = 1048576; /* (2^20) must be a power of two */
const unsigned int SOR_SIZE = 1000;    /*  NxN grid  */
const unsigned int SPARSE_SIZE_M = 100000;
const unsigned int SPARSE_SIZE_nz = 1000000;
const unsigned int LU_SIZE = 1000;

#endif
