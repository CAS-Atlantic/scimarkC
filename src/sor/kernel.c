/**
 *  SciMark2: A Java numerical benchmark measuring performance
 *  of computational kernels for FFTs, Monte Carlo simulation,
 *  sparse matrix computations, Jacobi SOR, and dense LU matrix
 *  factorizations.
 *
 *  Authors: Roldan Pozo (pozo@nist.gov) and Bruce Miller (bruce.miller@nist.gov)
 *
 *  Downloaded from: https://math.nist.gov/scimark2/download.html
 *
 *  Modified by: Aaron Graham (aaron.graham@unb.ca, aarongraham9@gmail.com) and
 *                Jean-Philippe Legault (jlegault@unb.ca, jeanphilippe.legault@gmail.com)
 *                 for the Centre for Advanced Studies - Atlantic (CAS-Atlantic) at the
 *                  Univerity of New Brunswick in Fredericton, New Brunswick, Canada
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "Random.h"
#include "Stopwatch.h"
#include "array.h"
#include "constants.h"

#include "kernel.h"

static void SOR_execute(int M,
                        int N,
                        double omega,
                        double** G,
                        int num_iterations) {
  double omega_over_four = omega * 0.25;
  double one_minus_omega = 1.0 - omega;

  /* update interior points */

  int Mm1 = M - 1;
  int Nm1 = N - 1;
  int p;
  int i;
  int j;
  double* Gi;
  double* Gim1;
  double* Gip1;

  for (p = 0; p < num_iterations; p++) {
    for (i = 1; i < Mm1; i++) {
      Gi = G[i];
      Gim1 = G[i - 1];
      Gip1 = G[i + 1];
      for (j = 1; j < Nm1; j++)
        Gi[j] = omega_over_four * (Gim1[j] + Gip1[j] + Gi[j - 1] + Gi[j + 1]) +
                one_minus_omega * Gi[j];
    }
  }
}

static double kernel_measureSOR(unsigned int N, int itter, Random R) {
  double** G = RandomMatrix(N, N, R);
  double result = 0.0;
  unsigned int i = 0, j = 0;

  Stopwatch Q = new_Stopwatch();
  Stopwatch_start(Q);
  SOR_execute(N, N, 1.25, G, itter);
  Stopwatch_stop(Q);
  result = Stopwatch_read(Q);
  Array2D_double_delete(N, N, G);
  Stopwatch_delete(Q);

  return result;
}

void kernel_measure(int itter, Random R) {
  printf("benchmark: SOR \nitteration: %d\ntime: %0.2fms\n", itter,
         kernel_measureSOR(SOR_SIZE, itter, R));
  fflush(stdout);
}