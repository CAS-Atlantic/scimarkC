#ifndef KERNEL_H
#define KERNEL_H

#include "Random.h"

double kernel_measureFFT(unsigned int FFT_size, int itter, Random R);

double kernel_measureSOR(unsigned int SOR_size, int itter, Random R);

double kernel_measureMonteCarlo(int itter, Random R);

double kernel_measureSparseMatMult(unsigned int Sparse_size_N,
                                   unsigned int Sparse_size_nz,
                                   int itter,
                                   Random R);

double kernel_measureLU(unsigned int LU_size, int itter, Random R);

#endif
