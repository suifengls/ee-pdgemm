/*Matrix-matrix Multipication with Pipeline Broadcast*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <mpi.h>
#include <unistd.h>
/* User defined headers */
#include "/home/lchen/myheader/pipeline.h"
#include "/home/lchen/myheader/map.h"

extern void pdgemm(char *TRANSA, char *TRANSB,
            int *M, int *N, int *K,
            double *ALPHA,
            double *A, int *IA, int *JA, int *DESCA,
            double *B, int *IB, int *JB, int *DESCB,
            double *BETA,
            double *C, int *IC, int *JC, int *DESCC);

extern void Cblacs_gridinfo( int context, int *nproc, int *npcol, int *myrow, int *mycol);

/* Local size of matrix, global size is 8*MATRIX_SIZE */
#define PIPE_SIZE 2000
/* DVFS status */
#define DVFS_HIGH 1
#define DVFS_LOW 0

// Descriptor
#define DTYPE 0
#define DCTXT 1
#define DM 2
#define DN 3
#define DMB 4
#define DNB 5
#define DRSRC 6
#define DCSRC 7
#define DLLD 8

/* processors in one node tiling is R or C major */
const int DVFS_ENABLE = 0;
/* Pipeline (1) or MPI(0) ? */
const int PIPE = 1;
/* timer inside the mmult? */
const int TIMER = 0;


/* number of processor rows and columns */
int rows;
int columns;
int ctxt;

/* External functions */
/* Dgemm function from ATLAS for local matrix-matrix multiplication */
extern void dgemm_(char*, char*, int*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);
/* Pipeline broadcast */
extern void MY_Bcast(double*, int, MPI_Datatype, int, MPI_Comm, int, int, int );
/* processor mapping for DVFS */
extern void mapping(int, int);

void pdgemm(char *TRANSA, char *TRANSB,
            int *M, int *N, int *K,
            double *ALPHA,
            double *a, int *IA, int *JA, int *DESCA,
            double *b, int *IB, int *JB, int *DESCB,
            double *BETA,
            double *c, int *IC, int *JC, int *DESCC)
{
   int m, n, k;
   int i, j, round;
   int r_rank, c_rank, g_rank, g_size;
   double *tempA, *tempB;
   double t_bc, t_calc, t_tot, t_bar;
   double tt_bc, tt_calc, tt_tot, tt_bar;
   double sum, g_sum;
   int roots;
   int pk;
   int nprow, npcol, myrow, mycol;
   MPI_Comm my_row_comm, my_column_comm;

   /* Get the rank of the global communicator */
   MPI_Comm_rank(MPI_COMM_WORLD, &g_rank);
   MPI_Comm_size(MPI_COMM_WORLD, &g_size);

   Cblacs_gridinfo( (ctxt = DESCC[DCTXT]), &nprow, &npcol, &myrow, &mycol);

   rows = nprow;
   columns = npcol;

   pk = PIPE_SIZE;

      m = DESCA[DM]/nprow;
      n = DESCB[DN]/nprow;
      k = DESCA[DN]/npcol; // local m n k size

   MPI_Barrier(MPI_COMM_WORLD);

   /* Create my row and column communicators */   
   MPI_Comm_split(MPI_COMM_WORLD, mycol, myrow, &my_column_comm);
   MPI_Comm_split(MPI_COMM_WORLD, myrow, mycol, &my_row_comm);

   /* Get the rank of the local row communicator */
   MPI_Comm_rank(my_row_comm, &r_rank);
   /* Get the rank of the local column communicator */
   MPI_Comm_rank(my_column_comm, &c_rank);
   
   tempA = (double *) malloc (m*k * sizeof(double));
   tempB = (double *) malloc (k*n * sizeof(double));

   //const enum CBLAS_ORDER Order=CblasRowMajor;
   //const enum CBLAS_TRANSPOSE TA=CblasNoTrans;
   //const enum CBLAS_TRANSPOSE TB=CblasNoTrans;

   sum = g_sum = 0.0;
   tt_bc = tt_calc = tt_tot = t_bar = 0.0;
   tt_tot = tt_bc = tt_calc = tt_bar = 0.0;

   if(TIMER)
      t_tot = MPI_Wtime();

      for(round=0; round<nprow; round++)
      {
         roots = round;

         if(r_rank == roots)
         {
              for(i=0; i<m; i++)
	        for(j=0; j<k; j++)
		{
		  tempA[i*k+j] = a[i*k+j];
		}
	 }

         if(c_rank == roots)
         {
              for(i=0; i<k; i++)
	        for(j=0; j<n; j++)
		{
		  tempB[i*n+j] = b[i*n+j];
		}
         }


         if(TIMER)
            t_bc = MPI_Wtime();

         /* Broadcast to right */
         if(PIPE)
	 {
              MY_Bcast(tempA, m*k,  MPI_DOUBLE, roots, my_row_comm, r_rank, columns, pk);
	 }
         else
         {
              MPI_Bcast(tempA, k*m, MPI_DOUBLE, roots, my_row_comm);
	 }
	 
         /* Broadcast below */
         
         if(PIPE)
	 {
              MY_Bcast(tempB, k*n, MPI_DOUBLE, roots, my_column_comm, c_rank, rows, pk);
	 }
         else
	 {
               MPI_Bcast(tempB, n*k, MPI_DOUBLE, roots, my_column_comm);
         }
	 
         if(TIMER)
         {
            t_bc = MPI_Wtime() - t_bc;
            tt_bc += t_bc;
            t_calc = MPI_Wtime();
         }

         if(DVFS_ENABLE) 
            mapping(g_rank%8, DVFS_HIGH);

         /* Do the multiplication */
         dgemm_(TRANSA, TRANSB, &m, &n, &k, ALPHA, tempA, &m, tempB, &k, BETA, c, &m);

         if(DVFS_ENABLE) 
            mapping(g_rank%8, DVFS_LOW);

         if(TIMER)
         { 
            t_calc = MPI_Wtime() - t_calc;
            tt_calc += t_calc;
         }
      }
         
   if(TIMER)
   {
      t_tot = MPI_Wtime() - t_tot;
      tt_tot += t_tot;
   
      MPI_Reduce(&tt_tot, &t_tot, 1, MPI_DOUBLE, MPI_MAX,  0, MPI_COMM_WORLD);
      MPI_Reduce(&tt_bc, &t_bc, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
      MPI_Reduce(&tt_calc, &t_calc, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

   
      if(g_rank == 0)
         printf("tot = %10.6f, bc = %10.6f, calc = %10.6f\n", t_tot, t_bc, t_calc);
   }

   MPI_Barrier(MPI_COMM_WORLD);

   free(tempA); free(tempB);
   MPI_Comm_free(&my_row_comm); 
   MPI_Comm_free(&my_column_comm);
}

