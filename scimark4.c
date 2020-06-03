#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Random.h"
#include "constants.h"
#include "kernel.h"

static void print_banner() {
  printf(
      "**                                                              **\n");
  printf(
      "** SciMark4 Numeric Benchmark, see http://math.nist.gov/scimark **\n");
  printf(
      "** for details. (Results can be submitted to pozo@nist.gov)     **\n");
  printf(
      "**                                                              **\n");
}

int main(int argc, char* argv[]) {
  /* default to the (small) cache-contained version */

  int itter = RESOLUTION_DEFAULT;
  Random R = new_Random_seed(RANDOM_SEED);

  if (argc > 1) {
    if (strcmp(argv[1], "-help") == 0 || strcmp(argv[1], "-h") == 0) {
      fprintf(stderr, "Usage: [itter]\n");
      exit(0);
    } else {
      itter = atoi(argv[1]);
    }
  }

  print_banner();
  printf("Using %d itter per kernel.", itter);
  printf("\n\n");

  /* print out results  */
  printf("FFT reps:              %0.2f\n",
         kernel_measureFFT(FFT_SIZE, itter, R));
  fflush(stdout);

  printf("SOR reps:              %0.2f\n",
         kernel_measureSOR(SOR_SIZE, itter, R));
  fflush(stdout);

  printf("Montel Carlo reps:     %0.2f\n",
         kernel_measureMonteCarlo(itter * MONTECARLO_ITTER, R));
  fflush(stdout);

  printf("Sparse MatMult repss:  %0.2f\n",
         kernel_measureSparseMatMult(SPARSE_SIZE_M, SPARSE_SIZE_nz, itter, R));
  fflush(stdout);

  printf("LU reps:               %0.2f\n", kernel_measureLU(LU_SIZE, itter, R));
  fflush(stdout);

  printf("\n");

  Random_delete(R);

  return 0;
}
