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

#ifndef STOPWATCH_H
#define STOPWATCH_H

#include <time.h>

typedef struct {
  int running; /* boolean */
  double last_time;
  double total;

} * Stopwatch, Stopwatch_struct;

double seconds();

void Stopwtach_reset(Stopwatch Q);

Stopwatch new_Stopwatch(void);
void Stopwatch_delete(Stopwatch S);
void Stopwatch_start(Stopwatch Q);
void Stopwatch_resume(Stopwatch Q);
void Stopwatch_stop(Stopwatch Q);
double Stopwatch_read(Stopwatch Q);

#endif
