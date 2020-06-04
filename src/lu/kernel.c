#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "Random.h"
#include "Stopwatch.h"
#include "array.h"
#include "constants.h"

#include "kernel.h"

static int LU_factor(int M, int N, double** A, int* pivot) {
  int minMN = M < N ? M : N;
  int j = 0;

  for (j = 0; j < minMN; j++) {
    /* find pivot in column j and  test for singularity. */

    int jp = j;
    int i;

    double t = fabs(A[j][j]);
    for (i = j + 1; i < M; i++) {
      double ab = fabs(A[i][j]);
      if (ab > t) {
        jp = i;
        t = ab;
      }
    }

    pivot[j] = jp;

    /* jp now has the index of maximum element  */
    /* of column j, below the diagonal          */

    if (A[jp][j] == 0)
      return 1; /* factorization failed because of zero pivot */

    if (jp != j) {
      /* swap rows j and jp */
      double* tA = A[j];
      A[j] = A[jp];
      A[jp] = tA;
    }

    if (j < M - 1) /* compute elements j+1:M of jth column  */
    {
      /* note A(j,j), was A(jp,p) previously which was */
      /* guarranteed not to be zero (Label #1)         */

      double recp = 1.0 / A[j][j];
      int k;
      for (k = j + 1; k < M; k++)
        A[k][j] *= recp;
    }

    if (j < minMN - 1) {
      /* rank-1 update to trailing submatrix:   E = E - x*y; */
      /* E is the region A(j+1:M, j+1:N) */
      /* x is the column vector A(j+1:M,j) */
      /* y is row vector A(j,j+1:N)        */

      int ii;
      for (ii = j + 1; ii < M; ii++) {
        double* Aii = A[ii];
        double* Aj = A[j];
        double AiiJ = Aii[j];
        int jj;
        for (jj = j + 1; jj < N; jj++)
          Aii[jj] -= AiiJ * Aj[jj];
      }
    }
  }

  return 0;
}

static double kernel_measureLU(unsigned int N, int itter, Random R) {
  double** A = NULL;
  double** lu = NULL;
  int* pivot = NULL;

  Stopwatch Q = new_Stopwatch();
  double result = 0.0;
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

  /* make sure A is diagonally dominant, to astatic void singularity */
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
    /* to astatic void being optimized out. */
    /*   lu_center = max( A_center, old_lu_center) */

    lu[N2][N2] = (A[N2][N2] > lu_center ? A[N2][N2] : lu_center);

    LU_factor(N, N, lu, pivot);
  }
  Stopwatch_stop(Q);

  /* approx Mflops */
  result = Stopwatch_read(Q);

  Stopwatch_delete(Q);
  free(pivot);

  Array2D_double_delete(N, N, lu);
  Array2D_double_delete(N, N, A);

  return result;
}

void kernel_measure(int itter, Random R) {
  printf("benchmark: LU \nitteration: %d\ntime: %0.2fms\n", itter,
         kernel_measureLU(LU_SIZE, itter, R));
  fflush(stdout);
}