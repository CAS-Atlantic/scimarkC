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
#include <stdio.h>
#include "array.h"

#ifndef NULL
#define NULL 0
#endif


double** new_Array2D_double(unsigned int M, unsigned int N)
{
    unsigned int i=0;
    int failed = 0;

    double **A = (double**) malloc(sizeof(double*)*M);
    if (A == NULL)
        return NULL;

    for (i=0; i<M; i++)
    {
        A[i] = (double*) malloc(N * sizeof(double));
        if (A[i] == NULL)
        {
            failed = 1;
            break;
        }
    }

    /* if we didn't successfully allocate all rows of A      */
    /* clean up any allocated memory (i.e. go back and free  */
    /* previous rows) and return NULL                        */

    if (failed)
    {
        i--;
        for (; i<=0; i--)
            free(A[i]);
        free(A);
        return NULL;
    }
    else
        return A;
}
void Array2D_double_delete(unsigned int M, unsigned int N, double **A)
{
    unsigned int i;
    if (A == NULL) return;

    for (i=0; i<M; i++)
        free(A[i]);

    free(A);
}


  void Array2D_double_copy(unsigned int M, unsigned int N, double **B, 
          double **A)
  {
        unsigned int remainder = N & 3;       /* N mod 4; */
        unsigned int i=0;
        unsigned int j=0;

        for (i=0; i<M; i++)
        {
            double *Bi = B[i];
            double *Ai = A[i];
            for (j=0; j<remainder; j++)
                Bi[j] = Ai[j];
            for (j=remainder; j<N; j+=4)
            {
                Bi[j] = Ai[j];
                Bi[j+1] = Ai[j+1];
                Bi[j+2] = Ai[j+2];
                Bi[j+3] = Ai[j+3];
            }
        }
  }
