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

static void SparseCompRow_matmult(int M,
                                  double* y,
                                  double* val,
                                  int* row,
                                  int* col,
                                  double* x,
                                  int NUM_ITERATIONS) {
  int reps;
  int r;
  int i;

  for (reps = 0; reps < NUM_ITERATIONS; reps++) {
    for (r = 0; r < M; r++) {
      double sum = 0.0;
      int rowR = row[r];
      int rowRp1 = row[r + 1];
      for (i = rowR; i < rowRp1; i++)
        sum += x[col[i]] * val[i];
      y[r] += sum;
    }
  }
}

static double kernel_measureSparseMatMult(unsigned int N,
                                          unsigned int nz,
                                          int itter,
                                          Random R) {
  /* initialize vector multipliers and storage for result */
  /* y = A*y;  */

  double* x = RandomVector(N, R);
  double* y = (double*)malloc(sizeof(double) * N);

  double result = 0.0;

#if 0
        // initialize square sparse matrix
        //
        // for this test, we create a sparse matrix with M/nz nonzeros
        // per row, with spaced-out evenly between the begining of the
        // row to the main diagonal.  Thus, the resulting pattern looks
        // like
        //             +-----------------+
        //             +*                +
        //             +***              +
        //             +* * *            +
        //             +** *  *          +
        //             +**  *   *        +
        //             +* *   *   *      +
        //             +*  *   *    *    +
        //             +*   *    *    *  + 
        //             +-----------------+
        //
        // (as best reproducible with integer artihmetic)
        // Note that the first nr rows will have elements past
        // the diagonal.
#endif

  int nr = nz / N;  /* average number of nonzeros per row  */
  int anz = nr * N; /* _actual_ number of nonzeros         */

  double* val = RandomVector(anz, R);
  int* col = (int*)malloc(sizeof(int) * nz);
  int* row = (int*)malloc(sizeof(int) * (N + 1));
  int r = 0;
  int i = 0;

  Stopwatch Q = new_Stopwatch();

  for (i = 0; i < N; i++)
    y[i] = 0.0;

  row[0] = 0;
  for (r = 0; r < N; r++) {
    /* initialize elements for row r */

    int rowr = row[r];
    int step = r / nr;
    int i = 0;

    row[r + 1] = rowr + nr;
    if (step < 1)
      step = 1; /* take at least unit steps */

    for (i = 0; i < nr; i++)
      col[rowr + i] = i * step;
  }

  Stopwatch_start(Q);
  SparseCompRow_matmult(N, y, val, row, col, x, itter);
  Stopwatch_stop(Q);

  /* approx Mflops */
  result = Stopwatch_read(Q);

  Stopwatch_delete(Q);
  free(row);
  free(col);
  free(val);
  free(y);
  free(x);

  return result;
}

void kernel_measure(int itter, Random R) {
  printf("benchmark: Sparse MatMult\nitteration: %d\ntime: %0.2fms\n", itter,
         kernel_measureSparseMatMult(SPARSE_SIZE_M, SPARSE_SIZE_nz, itter, R));
  fflush(stdout);
}