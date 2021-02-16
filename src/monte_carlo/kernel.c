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

static double MonteCarlo_integrate(int Num_samples, Random R);
static double MonteCarlo_num_flops(int Num_samples);

/**
 Estimate Pi by approximating the area of a circle.

 How: generate N random numbers in the unit square, (0,0) to (1,1)
 and see how are within a radius of 1 or less, i.e.
 <pre>

 sqrt(x^2 + y^2) < r

 </pre>
  since the radius is 1.0, we can square both sides
  and astatic void a sqrt() computation:
  <pre>

    x^2 + y^2 <= 1.0

  </pre>
  this area under the curve is (Pi * r^2)/ 4.0,
  and the area of the unit of square is 1.0,
  so Pi can be approximated by
  <pre>
                # points with x^2+y^2 < 1
     Pi =~      --------------------------  * 4.0
                     total # points

  </pre>

*/

static double MonteCarlo_num_flops(int Num_samples) {
  /* 3 flops in x^2+y^2 and 1 flop in random routine */

  return ((double)Num_samples) * 4.0;
}

static double MonteCarlo_integrate(int Num_samples, Random R) {
  int under_curve = 0;
  int count;

  for (count = 0; count < Num_samples; count++) {
    double x = Random_nextDouble(R);
    double y = Random_nextDouble(R);

    if (x * x + y * y <= 1.0)
      under_curve++;
  }

  return ((double)under_curve / Num_samples) * 4.0;
}

static double kernel_measureMonteCarlo(int itter, Random R) {
  double result = 0.0;
  Stopwatch Q = new_Stopwatch();

  Stopwatch_start(Q);
  MonteCarlo_integrate(itter, R);
  Stopwatch_stop(Q);

  /* approx Mflops */
  result = Stopwatch_read(Q);
  Stopwatch_delete(Q);
  return result;
}

void kernel_measure(int itter, Random R) {
  printf("benchmark: Monte Carlo\nitteration: %d\ntime: %0.2fms\n", itter,
         kernel_measureMonteCarlo(itter * MONTECARLO_ITTER, R));
}