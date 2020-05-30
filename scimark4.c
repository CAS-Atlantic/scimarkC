#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Random.h"
#include "constants.h"
#include "kernel.h"

void print_banner(void);

unsigned int ilog2(unsigned int N) {
  unsigned int exp = 0;
  while ((N = (N >> 1)))
    exp++;

  return exp;
}

int main(int argc, char* argv[]) {
  /* default to the (small) cache-contained version */

  int itter = RESOLUTION_DEFAULT;

  int FFT_size = FFT_SIZE;
  int SOR_size = SOR_SIZE;
  int Sparse_size_M = SPARSE_SIZE_M;
  int Sparse_size_nz = SPARSE_SIZE_nz;
  int LU_size = LU_SIZE;

  /* run the benchmark */
  double res[6] = {0.0};
  double sum[6] = {0.0}; /* checksum */
  double estimate_time[6] = {0.0};
  Random R = new_Random_seed(RANDOM_SEED);

  if (argc > 1) {
    int current_arg = 1;

    if (strcmp(argv[1], "-help") == 0 || strcmp(argv[1], "-h") == 0) {
      fprintf(stderr, "Usage: [-large | -huge #MB] [minimum_time]\n");
      exit(0);
    }

    if (strcmp(argv[1], "-large") == 0) {
      FFT_size = LG_FFT_SIZE;
      SOR_size = LG_SOR_SIZE;
      Sparse_size_M = LG_SPARSE_SIZE_M;
      Sparse_size_nz = LG_SPARSE_SIZE_nz;
      LU_size = LG_LU_SIZE;

      current_arg++;
    }

    if (current_arg < argc) {
      if (atof(argv[current_arg]) > 0.0)
        itter = atoi(argv[current_arg]);
    }
  }

  print_banner();
  printf("Using %d itter per kenel.", itter);
  printf("\n\n");

  /* print out results  */

  kernel_measureFFT(FFT_size, itter, R, &res[1], &sum[1], &estimate_time[1]);
  printf("FFT             Mflops: %8.2f    (N=%d) \n", res[1], FFT_size);

  kernel_measureSOR(SOR_size, itter, R, &res[2], &sum[2], &estimate_time[2]);
  printf("SOR             Mflops: %8.2f    (%d x %d) \n", res[2], SOR_size,
         SOR_size);

  kernel_measureMonteCarlo(itter, R, &res[3], &sum[3], &estimate_time[3]);
  printf("MonteCarlo:     Mflops: %8.2f  \n", res[3]);

  kernel_measureSparseMatMult(Sparse_size_M, Sparse_size_nz, itter * 1000, R,
                              &res[4], &sum[4], &estimate_time[4]);
  printf("Sparse matmult  Mflops: %8.2f    (N=%d, nz=%d)  \n", res[4],
         Sparse_size_M, Sparse_size_nz);

  kernel_measureLU(LU_size, itter, R, &res[5], &sum[5], &estimate_time[5]);
  printf("LU              Mflops: %8.2f    (M=%d, N=%d) \n", res[5], LU_size,
         LU_size);

  res[0] = (res[1] + res[2] + res[3] + res[4] + res[5]) / 5;
  sum[0] = (sum[1] + sum[2] + sum[3] + sum[4] + sum[5]) / 5;

  printf("\n");
  printf("************************************\n");
  printf("Composite Score:       %8.2f\n", res[0]);
  printf("************************************\n");
  printf("\n");

  printf("FFT reps:              %0.2f\n", estimate_time[1]);
  printf("SOR reps:              %0.2f\n", estimate_time[2]);
  printf("Montel Carlo reps:     %0.2f\n", estimate_time[3]);
  printf("Sparse MatMult repss:  %0.2f\n", estimate_time[4]);
  printf("LU reps:               %0.2f\n", estimate_time[5]);
  printf("\n");
  printf("checksum:              %20.16e\n", sum[0]);

  Random_delete(R);

  return 0;
}

void print_banner() {
  printf(
      "**                                                              **\n");
  printf(
      "** SciMark4 Numeric Benchmark, see http://math.nist.gov/scimark **\n");
  printf(
      "** for details. (Results can be submitted to pozo@nist.gov)     **\n");
  printf(
      "**                                                              **\n");
}
