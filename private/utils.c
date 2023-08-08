#include <math.h>
#include <string.h>
#include "mex.h"
#include "utils.h"

/* y += A x. A is a dense col-major (as always) matrix. */
void matvec(double *A, int m, int n, const double *x, double *y)
{
  int i, j, os;

  os = 0;
  for(j = 0; j < n; j++) {
    for(i = 0; i < m; i++)
      y[i] += A[os+i]*x[j];
    os += m;
  }
}

