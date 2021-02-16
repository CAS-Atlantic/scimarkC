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

/*-----------------------------------------------------------------------*/

static int int_log2(int n);
static void FFT_transform(int N, double* data);
static void FFT_inverse(int N, double* data);
static void FFT_bitreverse(int N, double* data);

static int int_log2(int n) {
  int k = 1;
  int log = 0;
  for (/*k=1*/; k < n; k *= 2, log++)
    ;
  if (n != (1 << log)) {
    printf("FFT: Data length is not a power of 2!: %d ", n);
    exit(1);
  }
  return log;
}

static void FFT_transform_internal(int N, double* data, int direction) {
  int n = N / 2;
  int bit = 0;
  int logn;
  int dual = 1;

  if (n == 1)
    return; /* Identity operation! */
  logn = int_log2(n);

  if (N == 0)
    return;

  /* bit reverse the input data for decimation in time algorithm */
  FFT_bitreverse(N, data);

  /* apply fft recursion */
  /* this loop executed int_log2(N) times */
  for (bit = 0; bit < logn; bit++, dual *= 2) {
    double w_real = 1.0;
    double w_imag = 0.0;
    int a;
    int b;

    double theta = 2.0 * direction * PI / (2.0 * (double)dual);
    double s = sin(theta);
    double t = sin(theta / 2.0);
    double s2 = 2.0 * t * t;

    for (a = 0, b = 0; b < n; b += 2 * dual) {
      int i = 2 * b;
      int j = 2 * (b + dual);

      double wd_real = data[j];
      double wd_imag = data[j + 1];

      data[j] = data[i] - wd_real;
      data[j + 1] = data[i + 1] - wd_imag;
      data[i] += wd_real;
      data[i + 1] += wd_imag;
    }

    /* a = 1 .. (dual-1) */
    for (a = 1; a < dual; a++) {
      /* trignometric recurrence for w-> exp(i theta) w */
      {
        double tmp_real = w_real - s * w_imag - s2 * w_real;
        double tmp_imag = w_imag + s * w_real - s2 * w_imag;
        w_real = tmp_real;
        w_imag = tmp_imag;
      }
      for (b = 0; b < n; b += 2 * dual) {
        int i = 2 * (b + a);
        int j = 2 * (b + a + dual);

        double z1_real = data[j];
        double z1_imag = data[j + 1];

        double wd_real = w_real * z1_real - w_imag * z1_imag;
        double wd_imag = w_real * z1_imag + w_imag * z1_real;

        data[j] = data[i] - wd_real;
        data[j + 1] = data[i + 1] - wd_imag;
        data[i] += wd_real;
        data[i + 1] += wd_imag;
      }
    }
  }
}

static void FFT_bitreverse(int N, double* data) {
  /* This is the Goldrader bit-reversal algorithm */
  int n = N / 2;
  int nm1 = n - 1;
  int i = 0;
  int j = 0;
  for (; i < nm1; i++) {
    /*int ii = 2*i; */
    int ii = i << 1;

    /*int jj = 2*j; */
    int jj = j << 1;

    /* int k = n / 2 ; */
    int k = n >> 1;

    if (i < j) {
      double tmp_real = data[ii];
      double tmp_imag = data[ii + 1];
      data[ii] = data[jj];
      data[ii + 1] = data[jj + 1];
      data[jj] = tmp_real;
      data[jj + 1] = tmp_imag;
    }

    while (k <= j) {
      /*j = j - k ; */
      j -= k;

      /*k = k / 2 ;  */
      k >>= 1;
    }
    j += k;
  }
}

static void FFT_transform(int N, double* data) {
  FFT_transform_internal(N, data, -1);
}

static void FFT_inverse(int N, double* data) {
  int n = N / 2;
  double norm = 0.0;
  int i = 0;
  FFT_transform_internal(N, data, +1);

  /* Normalize */

  norm = 1 / ((double)n);
  for (i = 0; i < N; i++)
    data[i] *= norm;
}

static double kernel_measureFFT(unsigned int N, int itter, Random R) {
  /* initialize FFT data as complex (N real/img pairs) */

  int twiceN = 2 * N;
  double* x = RandomVector(twiceN, R);
  Stopwatch Q = new_Stopwatch();
  int i = 0;
  double result = 0.0;

  Stopwatch_start(Q);
  for (i = 0; i < itter; i++) {
    FFT_transform(twiceN, x); /* forward transform */
    FFT_inverse(twiceN, x);   /* backward transform */
  }
  Stopwatch_stop(Q);
  /* approx Mflops */

  result = Stopwatch_read(Q);
  Stopwatch_delete(Q);
  free(x);

  return result;
}

void kernel_measure(int itter, Random R) {
  printf("benchmark: FFT\nitteration: %d\ntime: %0.2fms\n", itter,
         kernel_measureFFT(FFT_SIZE, itter, R));
  fflush(stdout);
}