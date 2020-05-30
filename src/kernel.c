#include "kernel.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "FFT.h"
#include "LU.h"
#include "MonteCarlo.h"
#include "Random.h"
#include "SOR.h"
#include "SparseCompRow.h"
#include "Stopwatch.h"
#include "array.h"

void kernel_measureFFT(unsigned int N,
                       int itter,
                       Random R,
                       double* res,
                       double* sum_,
                       double* estimate_time) {
  /* initialize FFT data as complex (N real/img pairs) */

  int twiceN = 2 * N;
  double* x = RandomVector(twiceN, R);
  Stopwatch Q = new_Stopwatch();
  int i = 0;
  double result = 0.0;
  double sum = 0.0;

  Stopwatch_start(Q);
  for (i = 0; i < itter; i++) {
    FFT_transform(twiceN, x); /* forward transform */
    FFT_inverse(twiceN, x);   /* backward transform */
  }
  Stopwatch_stop(Q);
  /* approx Mflops */

  *estimate_time = Stopwatch_read(Q);
  result = FFT_num_flops(N) * itter / *estimate_time * 1.0e-6;
  Stopwatch_delete(Q);

  for (i = 0; i < twiceN; i++) {
    sum += x[i];
  }
  sum /= twiceN;

  free(x);

  *sum_ = sum;
  *res = result;
}

void kernel_measureSOR(unsigned int N,
                       int itter,
                       Random R,
                       double* res,
                       double* sum_,
                       double* estimate_time) {
  double** G = RandomMatrix(N, N, R);
  double result = 0.0;
  double sum = 0.0;
  unsigned int i = 0, j = 0;

  Stopwatch Q = new_Stopwatch();
  Stopwatch_start(Q);
  SOR_execute(N, N, 1.25, G, itter);
  Stopwatch_stop(Q);

  /* approx Mflops */

  *estimate_time = Stopwatch_read(Q);
  result = SOR_num_flops(N, N, itter) / *estimate_time * 1.0e-6;

  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      sum += G[i][j];
  sum /= (N * N);

  Array2D_double_delete(N, N, G);
  *res = result;
  *sum_ = sum;

  Stopwatch_delete(Q);
}

void kernel_measureMonteCarlo(int itter,
                              Random R,
                              double* res,
                              double* sum_,
                              double* estimate_time) {
  double result = 0.0;
  double sum = 0.0;
  Stopwatch Q = new_Stopwatch();

  Stopwatch_start(Q);
  sum += MonteCarlo_integrate(itter);
  Stopwatch_stop(Q);

  /* approx Mflops */
  *estimate_time = Stopwatch_read(Q);
  result = MonteCarlo_num_flops(itter) / *estimate_time * 1.0e-6;
  Stopwatch_delete(Q);
  *res = result;
  *sum_ = sum;
}

void kernel_measureSparseMatMult(unsigned int N,
                                 unsigned int nz,
                                 int itter,
                                 Random R,
                                 double* res,
                                 double* sum_,
                                 double* estimate_time) {
  /* initialize vector multipliers and storage for result */
  /* y = A*y;  */

  double* x = RandomVector(N, R);
  double* y = (double*)malloc(sizeof(double) * N);

  double result = 0.0;
  double sum = 0.0;

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
  *estimate_time = Stopwatch_read(Q);
  result = SparseCompRow_num_flops(N, nz, itter) / *estimate_time * 1.0e-6;

  for (i = 0; i < N; i++)
    sum += y[i];
  sum /= N;

  Stopwatch_delete(Q);
  free(row);
  free(col);
  free(val);
  free(y);
  free(x);

  *res = result;
  *sum_ = sum;
}

void kernel_measureLU(unsigned int N,
                      int itter,
                      Random R,
                      double* res,
                      double* sum_,
                      double* estimate_time) {
  double** A = NULL;
  double** lu = NULL;
  int* pivot = NULL;

  Stopwatch Q = new_Stopwatch();
  double result = 0.0;
  double sum = 0.0;
  int i = 0, j = 0;
  int N2 = N / 2;

  if ((A = RandomMatrix(N, N, R)) == NULL)
    exit(1);
  if ((lu = new_Array2D_double(N, N)) == NULL)
    exit(1);
  if ((pivot = (int*)malloc(N * sizeof(int))) == NULL)
    exit(1);

  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      lu[i][j] = 0.0;

  /* make sure A is diagonally dominant, to avoid singularity */
  /* set diagonal to be 4 times the absolute value of its row sum */
  for (i = 0; i < N; i++) {
    double row_sum = 0.0;
    /* compute row sum of absoluate values  */
    for (j = 0; j < N; j++)
      row_sum += fabs(A[i][j]);
    A[i][i] = 4 * row_sum;
  }

  Stopwatch_start(Q);
  for (i = 0; i < itter; i++) {
    double lu_center = fabs(lu[N2][N2]);
    Array2D_double_copy(N, N, lu, A);

    /* add modification to A, based on previous LU */
    /* to avoid being optimized out. */
    /*   lu_center = max( A_center, old_lu_center) */

    lu[N2][N2] = (A[N2][N2] > lu_center ? A[N2][N2] : lu_center);

    LU_factor(N, N, lu, pivot);
  }
  Stopwatch_stop(Q);

  /* approx Mflops */
  *estimate_time = Stopwatch_read(Q);
  result = LU_num_flops(N) * itter / *estimate_time * 1.0e-6;

  Stopwatch_delete(Q);
  free(pivot);

  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      sum += lu[i][j];
  sum /= (N * N);

  Array2D_double_delete(N, N, lu);
  Array2D_double_delete(N, N, A);

  *res = result;
  *sum_ = sum;
}
