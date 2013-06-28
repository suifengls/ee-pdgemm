#ifdef MYPDGEMM_H
#define MYPDGEMM_H

#include "PBblacs.h"

extern void Cblacs_get(int context, int what, int*val);

extern void Cblacs_gridinit(int *context, char *order, int nproc_rows, int nproc_cols);

extern void Cblacs_pcoord(int context, int p, int *my_proc_row, int *my_proc_col);

extern void Cblacs_exit(int doneflag);

extern void Cblacs_gridinfo( int *context, int *nproc, int *npcol, int *myrow, int *mycol);

extern void Cblacs_gridmap(int *context, int *usermap, int, int ,int);

extern void pdgemm(char *TRANSA, char *TRANSB,
            int *M, int *N, int *K,
            double *ALPHA,
            double *A, int *IA, int *JA, int *DESCA,
            double *B, int *IB, int *JB, int *DESCB,
            double *BETA,
            double *C, int *IC, int *JC, int *DESCC);

#endif
