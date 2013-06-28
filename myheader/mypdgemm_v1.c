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

#define SHOW 12
/* Local size of matrix, global size is 8*MATRIX_SIZE */
#define PIPE_SIZE 2000
#define BLOCK_SIZE 32
/* number of processors in one node to tile */
#define PROC_NODE 8
#define DISP_SIZE 5
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
const int ROWTILING = 1;
const int DVFS_ENABLE = 0;
/* Pipeline (1) or MPI(0) ? */
const int PIPE = 1;
/* display own communicators? */
const int DISPLAY = 1;
/* timer inside the mmult? */
const int TIMER = 0;


/* number of processor rows and columns */
int rows = PROC_NODE;
int columns = PROC_NODE;
int ctxt;

/* External functions */
/* Dgemm function from ATLAS for local matrix-matrix multiplication */
extern void dgemm_(char*, char*, int*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);
/* Pipeline broadcast */
extern void MY_Bcast(double*, int, MPI_Datatype, int, MPI_Comm, int, int, int );
/* processor mapping for DVFS */
extern void mapping(int, int);

/* User defined functions */
/* Tile processor in each node */ 
void CreatTiling(int *);
/* Display tiled communicators */
void DisplayCommunicator(int, int, MPI_Comm, MPI_Comm);
/* Create own communicators */
void CreateReduceCommunicator(int *, int *, MPI_Comm *, MPI_Comm *);
/* Matrix-matrix multiplication */
void mmult(char *TRANSA, char *TRABSB, int M, int N, int K, double *APLHA, double *BETA,
           double *A, int IA, int JA, int *DESCA, 
           double *B, int IB, int JB, int *DESCB,
           double *c, int IC, int JC, int *DESCC,
           MPI_Comm my_row_comm, MPI_Comm my_column_comm, int bk);

void CreatTiling(int *position)
{
   int myrank, nprocs;
   int p, r, c;
   int nc, nr, ns;
   p = 0;
   ns = PROC_NODE; 
   
   if(ROWTILING)
   { 
      nc = (int)(sqrt(ns));
      nr = ns/nc;
      while(nr*nc != ns)
      {
         nc--;
         nr = ns/nc;
      }
   }
   else
   { 
      nr = (int)(sqrt(ns));
      nc = ns/nr;
      while(nr*nc != ns)
      {
         nr--;
         nc = ns/nr;
      }
   }

   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   /*How many Global rows and columns?*/
   /*How many rows and columns?*/
   MPI_Barrier(MPI_COMM_WORLD);
   if(myrank == 0)
   {
      for(r=0; r<rows; r++)
      {
         for(c=0; c<columns; c++)
         {
            if(c%nc == 0)
               p = (r%nr)*nc+(r/nr)*(columns*nr)+(c/nc)*ns;
            else 
              p++;
            position[r*columns+c] = p;
         }
      } 
   }
   MPI_Bcast(position, nprocs, MPI_INT, 0, MPI_COMM_WORLD);

} 


void CreateReduceCommunicator(int *position, int *coordinate, MPI_Comm *my_row_comm, MPI_Comm *my_column_comm)
{
   int myproc, nprocs;
   int my_row, my_column;

   MPI_Group group_ordered, group_world;
   MPI_Comm ordered_comm;

   /* Access and set amount of prcessors */
   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
   /* Get the group underlying MPI_COMM_WORLD */
   MPI_Comm_group(MPI_COMM_WORLD, &group_world);
   /* Create the new group */
   MPI_Group_incl(group_world, nprocs, position, &group_ordered);
   /* Create the new Communicator */
   MPI_Comm_create(MPI_COMM_WORLD, group_ordered, &ordered_comm);
   /* Get the rank of the new communicator */
   MPI_Comm_rank(ordered_comm, &myproc);

   /* Create the Row and Columns Communicators */ 
   /* Assign my row and columns number */
   my_row = myproc/columns;
   my_column = myproc%columns;
   /* Create the row communicator */
   MPI_Comm_split(ordered_comm, my_row, myproc, my_row_comm);
   /* Create the columns comunicator */
   MPI_Comm_split(ordered_comm, my_column, myproc, my_column_comm);

   /* Save the coordinate information */
   coordinate[0] = my_row;
   coordinate[1] = my_column;

   //printf("position = %d, row = %d, column = %d\n", position[myproc],coordinate[0], coordinate[1]);
           

   /* Remove the ordered communicator */
   MPI_Group_free(&group_ordered);
   MPI_Comm_free(&ordered_comm);
}

void DisplayCommunicator(int my_row, int my_column, MPI_Comm my_row_comm, MPI_Comm my_column_comm)
{
   int i, myproc, *rank=NULL;
   int local_rank, local_p;

   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Comm_rank(MPI_COMM_WORLD, &myproc);
   /* Access and set amount of prcessors in local communicator */
   MPI_Comm_size(my_row_comm, &local_p);
   /* Get the rank of the local communicator */
   MPI_Comm_rank(my_row_comm, &local_rank);
   /* Get the members of the local communicator */
   if(local_rank == 0)
      rank = (int *)malloc(local_p*sizeof(int));

   MPI_Gather(&myproc, 1, MPI_INT, rank, 1, MPI_INT, 0, my_row_comm);

   MPI_Barrier(MPI_COMM_WORLD);

   if(local_rank == 0)
   {
      /* Display the Communicator */
      printf("Row   : %2d {%2d", my_row, rank[0]);
      for(i=1; i<local_p; i++)
         printf(", %2d", rank[i]);

      printf("}\n");
   }

   MPI_Barrier(MPI_COMM_WORLD);
   /* Access and set amount of prcessors in local communicator */
   MPI_Comm_size(my_column_comm, &local_p);
   /* Get the rank of the local communicator */
   MPI_Comm_rank(my_column_comm, &local_rank);

   /* Get the members of the local communicator */
   if(local_rank == 0)
      rank = (int *)malloc(local_p*sizeof(int));

   MPI_Gather(&myproc, 1, MPI_INT, rank, 1, MPI_INT, 0, my_column_comm);

   MPI_Barrier(MPI_COMM_WORLD);

   if(local_rank == 0)
   {
      /* Display the Communicator */
      printf("Column: %2d {%2d", my_column, rank[0]);
      for(i=1; i<local_p; i++)
         printf(", %2d", rank[i]);

      printf("}\n");
   }
}

void pdgemm(char *TRANSA, char *TRANSB,
            int *M, int *N, int *K,
            double *ALPHA,
            double *a, int *IA, int *JA, int *DESCA,
            double *b, int *IB, int *JB, int *DESCB,
            double *BETA,
            double *c, int *IC, int *JC, int *DESCC)
{
  //   int my_row, my_column;
   int bk = BLOCK_SIZE;
   int n;
   int i, j, k, round;
   int r_rank, c_rank, g_rank, g_size;
   double *tempA, *tempB;
   double *tmpA, *tmpB;
   double t_bc, t_calc, t_tot, t_bar;
   double tt_bc, tt_calc, tt_tot, tt_bar;
   double sum, g_sum;
   int num_bk, loops, roots;
   int pk, last=0;
   int nprow, npcol, myrow, mycol;
   int *position, *coordinate;
   MPI_Comm my_row_comm, my_column_comm;

   /* Get the rank of the global communicator */
   MPI_Comm_rank(MPI_COMM_WORLD, &g_rank);
   MPI_Comm_size(MPI_COMM_WORLD, &g_size);

   Cblacs_gridinfo( (ctxt = DESCC[DCTXT]), &nprow, &npcol, &myrow, &mycol);

   //   printf("nprow = %d, npcol = %d, myrow = %d, mycol = %d\n", nprow, npcol, myrow, mycol);

   pk = PIPE_SIZE;
   if(nprow == npcol)
      n = DESCA[DN]/nprow; // local size
   else 
   { 
     printf("The grid of process should be square.\n");
     exit(-1);
   }

   coordinate = (int *)malloc(2*sizeof(int));
   position = (int *)malloc(g_size*sizeof(int));
   

   //MPI_Allgather(pos, 1, MPI_INT, position, g_size, MPI_INT, MPI_COMM_WORLD);
   for(i=0; i<g_size; i++)
     position[i]=i;   
   MPI_Barrier(MPI_COMM_WORLD);
   //  CreatTiling(position);
   //CreateReduceCommunicator(position, coordinate, &my_row_comm, &my_column_comm);
   //my_row = coordinate[0];
   //my_column = coordinate[1];

   /* Create my row and column communicators */   
   MPI_Comm_split(MPI_COMM_WORLD, mycol, myrow, &my_column_comm);
   MPI_Comm_split(MPI_COMM_WORLD, myrow, mycol, &my_row_comm);

   /* Get the rank of the local row communicator */
   //MPI_Comm_rank(ctxtr->scp->comm, &r_rank);
   MPI_Comm_rank(my_row_comm, &r_rank);
   /* Get the rank of the local column communicator */
   MPI_Comm_rank(my_column_comm, &c_rank);
   //MPI_Comm_rank(ctxtc->scp->comm, &c_rank);
   if(r_rank == 0 )
     printf("my c_rank = %d\n", c_rank);


   tempA = (double *) malloc (n*bk * sizeof(double));
   tempB = (double *) malloc (bk*n * sizeof(double));

   //const enum CBLAS_ORDER Order=CblasRowMajor;
   //const enum CBLAS_TRANSPOSE TA=CblasNoTrans;
   //const enum CBLAS_TRANSPOSE TB=CblasNoTrans;

   sum = g_sum = 0.0;
   tt_bc = tt_calc = tt_tot = t_bar = 0.0;
   tt_tot = tt_bc = tt_calc = tt_bar = 0.0;
   num_bk =(n/bk);
   last = n%bk;

   tmpA = (double *) malloc (n*last * sizeof(double));
   tmpB = (double *) malloc (last*n * sizeof(double));

   if(last==0)
      loops = rows*num_bk;
   else 
   {
      num_bk++;
      loops = rows*num_bk;
   }
   

   if(TIMER)
      t_tot = MPI_Wtime();

      for(round=0; round<loops; round++)
      {
        k = (round%num_bk)*bk;
        roots = round/num_bk;


	if(*TRANSA == 'N')
	{
         if(r_rank == roots)
         {
           if(last!=0 && (round+1)%num_bk==0) // the last block of this process
	   {
              for(i=0; i<last; i++)
   	        for(j=0; j<n; j++)
		{
		  tmpA[i*n+j] = a[(k+i)*n+j];
		}
	   }
           else 
	   {
              for(i=0; i<bk; i++)
	        for(j=0; j<n; j++)
		{
		  tempA[i*n+j] = a[(k+i)*n+j];
		}
           }
	 }

         if(c_rank == roots)
         {
           if(last!=0 && (round+1)%num_bk==0) // the last block of this process
	   {
              for(i=0; i<n; i++)
	        for(j=0; j<last; j++)
		{
		  tmpB[i*last+j] = b[i*n+(k+j)];
		}
	   }
           else
	   {
              for(i=0; i<n; i++)
	        for(j=0; j<bk; j++)
		{
		  tempB[i*bk+j] = b[i*n+(k+j)];
		}
	   }
         }
	}
        else if( *TRANSA == 'T')
	{
         if(r_rank == roots)
         {
           if(last!=0 && (round+1)%num_bk==0) // the last block of this process
	   {
              for(i=0; i<n; i++)
   	        for(j=0; j<last; j++)
		{
		   tmpA[i*last+j] = a[i*n+(k+j)];
       		}
	   }
	   else
	   {
              for(i=0; i<n; i++)
	        for(j=0; j<bk; j++)
		{
		   tempA[i*bk+j] = a[i*n+(k+j)];
       		}
	   }
         }
         if(c_rank == roots)
         {
           if(last!=0 && (round+1)%num_bk==0) // the last block of this process
	   {       
              for(i=0; i<last; i++)
	        for(j=0; j<n; j++)
		{
		  tmpB[i*n+j] = b[(k+i)*n+j];
       		}
	   }
	   else
	   {
             for(i=0; i<bk; i++)
	        for(j=0; j<n; j++)
		{
		  tempB[i*n+j] = b[(k+i)*n+j];
       		}
	   }
         }
	}
   
         if(TIMER)
            t_bc = MPI_Wtime();

         /* Broadcast to right */
         if(PIPE)
	 {
           if(last!=0 && (round+1)%num_bk==0) // the last block of this process
              MY_Bcast(tmpA, n*last,  MPI_DOUBLE, roots, my_row_comm, r_rank, columns, pk);
	   else
              MY_Bcast(tempA, n*bk,  MPI_DOUBLE, roots, my_row_comm, r_rank, columns, pk);
	 }
         else
         {
           if(last!=0 && (round+1)%num_bk==0) // the last block of this process
              MPI_Bcast(tmpA, last*n, MPI_DOUBLE, roots, my_row_comm);
	   else
              MPI_Bcast(tempA, bk*n, MPI_DOUBLE, roots, my_row_comm);
	 }
	 
         /* Broadcast below */
         
         if(PIPE)
	 {
           if(last!=0 && (round+1)%num_bk==0) // the last block of this process
	      MY_Bcast(tmpB, last*n, MPI_DOUBLE, roots, my_column_comm, c_rank, rows, pk);
	   else
              MY_Bcast(tempB, bk*n, MPI_DOUBLE, roots, my_column_comm, c_rank, rows, pk);
	 }
         else
	 {
            if(last!=0 && (round+1)%num_bk==0) // the last block of this process
               MPI_Bcast(tmpB, n*last, MPI_DOUBLE, roots, my_column_comm);
	    else
               MPI_Bcast(tempB, n*bk, MPI_DOUBLE, roots, my_column_comm);
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
         if(*TRANSA == 'N') 
         {
           if(last!=0 && (round+1)%num_bk==0) // the last block of this process
              dgemm_(TRANSA, TRANSB, &n, &n, &last, ALPHA, tmpA, &n, tmpB, &last, BETA, c, &n);
	   else
              dgemm_(TRANSA, TRANSB, &n, &n, &bk, ALPHA, tempA, &n, tempB, &bk, BETA, c, &n);
	 }
         else if(*TRANSA == 'T') 
         {
           if(last!=0 && (round+1)%num_bk==0) // the last block of this process              
              dgemm_(TRANSA, TRANSB, &n, &n, &last, ALPHA, tmpA, &last, tmpB, &n, BETA, c, &n);
	   else
              dgemm_(TRANSA, TRANSB, &n, &n, &bk, ALPHA, tempA, &bk, tempB, &n, BETA, c, &n);
	 }

	 //	          dgemm_(TRANSA, TRANSB, &n, &n, &bk, ALPHA, tempA, &n, tempB, &bk, BETA, c, &n);
	 //	 cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, bk, alpha, tempA, bk, tempB, n, beta, c, n);

         /* check the calculate results*/
     
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

   free(tmpA);free(tmpB);
   //   free(a); free(b); free(c);
   free(tempA); free(tempB);
}

