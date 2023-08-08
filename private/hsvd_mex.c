/* Works with hsvd.m. */

#include <math.h>
#include <string.h>
#include <omp.h>
#include "utils.h"

typedef struct block_t {
  // rows = r:r+m-1; cols = c:c+n-1
  int r, c, m, n;
  // rank if U,Vt are used
  int k;
  // [U (with s mult'ed in) and V^T] or [dense block B]
  double *U, *Vt, *B;
} Block;

/* Convert the matlab block struct to the mex one.
     n: r c U Vt B
 */
void Block_New (Block *b, const mxArray *ns, int i)
{
  int s;
  double *tmp;

  tmp = mxGetPr(mxGetField(ns, i, "r"));
  b->r = (int)tmp[0] - 1;
  b->m = (int)(tmp[1] - tmp[0] + 1);
  tmp = mxGetPr(mxGetField(ns, i, "c"));
  b->c = (int)tmp[0] - 1;
  b->n = (int)(tmp[1] - tmp[0] + 1);
  if (mxIsEmpty(mxGetField(ns, i, "B"))) {
    b->k = mxGetN(mxGetField(ns, i, "U"));
    s = b->m*b->k;
    b->U = malloc(s*sizeof(double));
    memcpy(b->U,  mxGetPr(mxGetField(ns, i, "U")),  s*sizeof(double));
    s = b->k*b->n;
    b->Vt = malloc(s*sizeof(double));
    memcpy(b->Vt, mxGetPr(mxGetField(ns, i, "Vt")), s*sizeof(double));
    b->B = NULL;
  }
  else {
    s = b->m*b->n;
    b->B = malloc(s*sizeof(double));
    memcpy(b->B,  mxGetPr(mxGetField(ns, i, "B")),  s*sizeof(double));
    b->U = b->Vt = NULL;
  }
}

void Block_Free (Block *b)
{
  if (b->U) {
    free(b->U);
    free(b->Vt);
  }
  if (b->B) free(b->B);
}

void Block_matvec (const Block *b, const double *x, double *y, double *work)
{
  if (b->B) {
    // y(r(1):r(2)) = B x(c(1):c(2))
    matvec(b->B,  b->m, b->n, x + b->c, y + b->r);
  } else {
    memset(work, 0, b->k*sizeof(double));
    // work = V^T x(c(1):c(2))
    matvec(b->Vt, b->k, b->n, x + b->c, work);
    // y(r(1):r(2)) = U work
    matvec(b->U,  b->m, b->k, work,     y + b->r);
  }
}

typedef struct hsvd_t {
  int nthreads, m, *nblks, *par;
  Block **blk;
  double *work;
} Hsvd;

/* Sits in memory between mex calls. Can only hold two, the max number we need
   in our simulations. */
static Hsvd gh[8];
bool goccupied[8] = {false, false, false, false, false, false, false, false};

/* Convert the matlab hsvd struct to the mex one.
     hsvd: par sz err ns
 */
void Hsvd_New (Hsvd *h, const mxArray *mh)
{
  int i, j;
  mxArray *ns, *nsc;
  double *tmp;

  ns = mxGetField(mh, 0, "ns");
  h->nthreads = mxGetNumberOfElements(ns);
  h->nblks = malloc(h->nthreads*sizeof(int));
  h->blk = malloc(h->nthreads*sizeof(Block*));
  for (i = 0; i < h->nthreads; i++) {
    nsc = mxGetCell(mxGetField(mh, 0, "ns"), i);
    h->nblks[i] = mxGetNumberOfElements(nsc);
    h->blk[i] = malloc(h->nblks[i]*sizeof(Block));
    for (j = 0; j < h->nblks[i]; j++)
      Block_New(h->blk[i] + j, nsc, j);
  }

  h->par = malloc(h->nthreads*sizeof(int));
  tmp = mxGetPr(mxGetField(mh, 0, "par"));
  for (i = 0; i < h->nthreads; i++)
    h->par[i] = (int)tmp[i] - 1;

  h->m = (int)mxGetPr(mxGetField(mh, 0, "sz"))[0];
  h->work = malloc(h->m*sizeof(double));
}

void Hsvd_Free (Hsvd *h)
{
  int i, j;

  for (i = 0; i < h->nthreads; i++) {
    for (j = 0; j < h->nblks[i]; j++)
      Block_Free(h->blk[i] + j);
    free(h->blk[i]);
  }
  free(h->blk);
  free(h->nblks);
  free(h->par);
  free(h->work);
}

/* matvec with hsvd struct */
void Hsvd_matvec (const Hsvd *h, const double *x, double *y)
{
  int tid, i, nthreads;

#pragma omp parallel private(tid,i) num_threads(h->nthreads)
  {
    tid = omp_get_thread_num();
    
    if (tid == 0)
      nthreads = omp_get_num_threads();
    
    for (i = 0; i < h->nblks[tid]; i++)
      Block_matvec(h->blk[tid] + i, x, y, h->work + h->par[tid]);
  }

  // Did we miss any blocks b/c the number of CPU threads was less than
  // expected?
  for (tid = nthreads; tid < h->nthreads; tid++)
    for (i = 0; i < h->nblks[tid]; i++)
      Block_matvec(h->blk[tid] + i, x, y, h->work + h->par[tid]);
}

void Cleanup (int id)
{
  if (goccupied[id]) Hsvd_Free(gh+id);
  goccupied[id] = false;
}

void Init (const mxArray *mh, int id)
{
  Cleanup(id);
  Hsvd_New(gh+id, mh);
  goccupied[id] = true;
}

/* This function works with hsvd.m and implements
     hsvd_mex(h,'init',id)
     hsvd_mex(0,'cleanup',h.id)
     y = hsvd_mex(h.id,x)
*/
void mexFunction (int nlhs, mxArray **plhs, int nrhs, const mxArray **prhs)
{
  int id;

  if (nrhs == 3) {
    char msg[2];

    mxGetString(prhs[1], msg, 2);
    if (msg[0] == 'i') {
      id = (int)mxGetScalar(prhs[2]);
      Init(prhs[0], id);
    } else {
      id = (int)mxGetScalar(prhs[0]);
      Cleanup(id);
    }
  } else {
    id = (int)mxGetScalar(prhs[0]);
    plhs[0] = mxCreateDoubleMatrix(gh[id].m, 1, mxREAL);
    Hsvd_matvec(gh+id, mxGetPr(prhs[1]), mxGetPr(plhs[0]));
  }
}
