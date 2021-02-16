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

#ifndef RANDOM_H
#define RANDOM_H

typedef struct {
  int m[17];
  int seed;
  int i;                     /* originally = 4 */
  int j;                     /* originally =  16 */
  int /*boolean*/ haveRange; /* = false; */
  double left;               /*= 0.0; */
  double right;              /* = 1.0; */
  double width;              /* = 1.0; */
} Random_struct, *Random;

Random new_Random_seed(int seed);
double Random_nextDouble(Random R);
void Random_delete(Random R);
double* RandomVector(unsigned int N, Random R);
double** RandomMatrix(unsigned int M, unsigned int N, Random R);

#endif