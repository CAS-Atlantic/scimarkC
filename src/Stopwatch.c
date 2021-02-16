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

#include <stdlib.h>
#include "Stopwatch.h"

double seconds()
{
        return ((double) clock()) / (double) CLOCKS_PER_SEC; 
}

void Stopwtach_reset(Stopwatch Q)
{
    Q->running = 0;         /* false */
    Q->last_time = 0.0;
    Q->total= 0.0;
}


Stopwatch new_Stopwatch(void)
{
    Stopwatch S = (Stopwatch) malloc(sizeof(Stopwatch_struct));
    if (S == NULL)
        return NULL;

    Stopwtach_reset(S);
    return S;
}

void Stopwatch_delete(Stopwatch S)
{
    if (S != NULL)
        free(S);
}


/* Start resets the timer to 0.0; use resume for continued total */

void Stopwatch_start(Stopwatch Q)
{
    if (! (Q->running)  )
    {
        Q->running = 1;  /* true */
        Q->total = 0.0;
        Q->last_time = seconds();
    }
}
   
/** 
    Resume timing, after stopping.  (Does not wipe out
        accumulated times.)

*/

void Stopwatch_resume(Stopwatch Q)
{
    if (!(Q->running))
    {
        Q-> last_time = seconds(); 
        Q->running = 1;  /*true*/
    }
}
   
void Stopwatch_stop(Stopwatch Q)  
{ 
    if (Q->running) 
    { 
        Q->total += seconds() - Q->last_time; 
        Q->running = 0;  /* false */
    }
}
  
 
double Stopwatch_read(Stopwatch Q)
{  
    
    if (Q->running) 
    {
        double t = seconds();
        Q->total += t - Q->last_time;
        Q->last_time = t;
    }
    return Q->total;
}
        
