/*Matrix-matrix Multipication with Pipeline Broadcast*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <mpi.h>
#include <unistd.h>

/* User defined headers */
#include "/home/lchen/pdgemm/pdgemm.h"
#include "/home/lchen/myheader/pipeline.h"
#include "/home/lchen/myheader/map.h"
#include "/home/lchen/myheader/mypdgemm.h"

extern void pdgemm(char *TRANSA, char *TRANSB,
            int *M, int *N, int *K,
            double *ALPHA,
            double *A, int *IA, int *JA, int *DESCA,
            double *B, int *IB, int *JB, int *DESCB,
            double *BETA,
            double *C, int *IC, int *JC, int *DESCC);

/* Local size of matrix, global size is 8*MATRIX_SIZE */
#define GLOBAL_M 26000
#define GLOBAL_N 26000
#define GLOBAL_K 26000
/* block size of pdgemm in each process */
#define EACH 65
/* number of processors in one node to tile */
#define PROC_NODE 8
#define MATRIX_SIZE ((GLOBAL_M)/PROC_NODE)
#define DESC_SIZE 9
#define DISP_SIZE 5
#define SHOW 61 // show id 0-63, 64 - not show
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

/* meansure energy? */
const int MEASURE = 1;
/* display hostname? */
const int HELLO = 0;
/* major 0 - not tiling, 1 - row major tiling, 2- column major tiling */
const int MAJOR = 1;

/* External functions */
/* processor mapping for DVFS */
extern void mapping(int, int);

/* Create random matrix */
void RndMatrix(double*, int, int, int);
/* Create zero matrix */
void Zeros(double*, int, int);
/* Copy from a to b */
void CpyMatrix(double *a, double *b, int r, int c);
/* The normal difference between a and b */
double diff_norm(double* a, double* b, int, int);

int main (int argc, char **argv)
{
   int local_mat_rows;
   int local_mat_cols;
   double *A_local;
   int A_descrip[DESC_SIZE];
   double *B_local;
   int B_descrip[DESC_SIZE];
   double *C_local;
   int C_descrip[DESC_SIZE];

   int nproc_rows;
   int nproc_cols;
   int m, n;
   int row_block_size;
   int col_block_size;
   int blacs_grid;
   int my_process_row;
   int my_process_col;

   int myproc, nprocs;
   char myname[MPI_MAX_PROCESSOR_NAME];
   double *a, *b, *c;
   MPI_Comm my_row_comm, my_column_comm;

   int M, N, K;
   //   int local_m, local_n, local_k;

   M = GLOBAL_M;
   N = GLOBAL_N;
   K = GLOBAL_K;  

   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
   MPI_Comm_rank(MPI_COMM_WORLD, &myproc);

   /* Ensure we have at least two processors */
   if (nprocs < 2) 
   {
       printf("Too few processors!\n");
       exit (1);
   }

   if(gethostname (myname, MPI_MAX_PROCESSOR_NAME) != 0 )
      printf("Error: gethostname failed!\n");
   else if(HELLO)
      printf("Hello from %2d of %2d on %s\n", myproc, nprocs, myname);

   /* Set to HIGH frequency */
   mapping(myproc%8, DVFS_HIGH);

   Get_input(nprocs, myproc, &n, &nproc_rows, &nproc_cols, &row_block_size, &col_block_size);
   m = n;

   Cblacs_get(0, 0, &blacs_grid);
   /* ROW MAJOR TILING */
   if(MAJOR==1)
   {
      int usermap[64]=  {0,  1,  8,  9,  16, 17, 24, 25,
                       2,  3,  10, 11, 18, 19, 26, 27,
                       4,  5,  12, 13, 20, 21, 28, 29,
                       6,  7,  14, 15, 22, 23, 30, 31,
                       32, 33, 40, 41, 48, 49, 56, 57,
                       34, 35, 42, 43, 50, 51, 58, 59,
                       36, 37, 44, 45, 52, 53, 60, 61,
                       38, 39, 46, 47, 54, 55, 62, 63};
      int ldumap=8, nprow=8, npcol=8;
      Cblacs_gridmap(&blacs_grid, usermap, ldumap, nprow, npcol);
   }
   else if (MAJOR==2)	
    /* COLUMN MAJOR TILING*/
   {
      int usermap[64]={0, 1, 2, 3, 8, 9, 10, 11,
                       4, 5, 6, 7, 12,13, 14, 15,
                       16, 17, 18, 19, 24, 25, 26, 27,
                       20, 21, 22, 23, 28, 29, 30, 31,
                       32, 33, 34, 35, 40, 41, 42, 43,
                       36, 37, 38, 39, 44, 45, 46, 47,
                       48, 49, 50, 51, 56, 57, 58, 59,
                       52, 53, 54, 55, 60, 61, 62, 63};
    
      int ldumap=8, nprow=8, npcol=8;
      Cblacs_gridmap(&blacs_grid, usermap, ldumap, nprow, npcol);
   }
   else if(MAJOR==0)
      Cblacs_gridinit(&blacs_grid, "R", nproc_rows, nproc_cols);
   
   Cblacs_pcoord(blacs_grid, myproc, &my_process_row, &my_process_col);

   local_mat_rows = Get_dimension(m, row_block_size, my_process_row, nproc_rows);
   local_mat_cols = Get_dimension(n, col_block_size, my_process_col, nproc_cols);


   a = (double *) malloc (local_mat_rows*local_mat_cols * sizeof(double));
   b = (double *) malloc (local_mat_rows*local_mat_cols * sizeof(double));
   c = (double *) malloc (local_mat_rows*local_mat_cols * sizeof(double));

   A_local = (double *) malloc (local_mat_rows*local_mat_cols * sizeof(double));
   B_local = (double *) malloc (local_mat_rows*local_mat_cols * sizeof(double));
   C_local = (double *) malloc (local_mat_rows*local_mat_cols * sizeof(double));


   Build_descrip(myproc, "A", A_descrip, m, n, row_block_size, col_block_size, blacs_grid, local_mat_rows);
   Build_descrip(myproc, "B", B_descrip, m, n, row_block_size, col_block_size, blacs_grid, local_mat_rows);
   Build_descrip(myproc, "C", C_descrip, m, n, row_block_size, col_block_size, blacs_grid, local_mat_rows);

   if(A_descrip[DM]%A_descrip[DMB]!=0)
   {
     printf("The block size should be divideable by the global size.\n");
     exit(-1);
   }

   if(myproc==SHOW)
     printf("\nA_descrip = [ %d, %d, %d, %d, %d, %d, %d, %d, %d]\n", 
     A_descrip[0], A_descrip[1], A_descrip[2], A_descrip[3], A_descrip[4], A_descrip[5], A_descrip[6], A_descrip[7], A_descrip[8]);

   RndMatrix(a, local_mat_rows, local_mat_cols, myproc); 
   CpyMatrix(a, A_local, local_mat_rows, local_mat_cols);
   RndMatrix(b, local_mat_rows, local_mat_cols, myproc); 
   CpyMatrix(b, B_local, local_mat_rows, local_mat_cols);
   Zeros(c, local_mat_rows, local_mat_cols); 
   Zeros(C_local, local_mat_rows, local_mat_cols); 
       
   int ij = 1;
   char tran = 'N';
   double alpha = 1.0, beta = 1.0;
   double exetime, exetime0;

   MPI_Barrier(MPI_COMM_WORLD);
  
   if(MEASURE && myproc==0)
   {
      system("/apps/power-bench/mclient -H 10.1.255.100 -d /tmp");
      system("/apps/power-bench/mclient -H 10.1.255.100 -l pdgemm.ptr");
      system("/apps/power-bench/mclient -H 10.1.255.100 -s pdgemm");
   }

 
   MPI_Barrier(MPI_COMM_WORLD);
   exetime = -MPI_Wtime();
   /* My pdgemm */
   pdgemm(&tran, &tran, &M, &N, &K, &alpha, a, &ij, &ij, A_descrip, b, &ij, &ij, B_descrip, &beta, c, &ij, &ij, C_descrip);    

   MPI_Barrier(MPI_COMM_WORLD);
   exetime += MPI_Wtime();

   MPI_Barrier(MPI_COMM_WORLD);
   exetime0 = -MPI_Wtime();
   /* ScaLAPACK pdgemm */
   Mat_vect_mult(m, n, A_local, A_descrip, B_local, B_descrip, C_local, C_descrip);

   MPI_Barrier(MPI_COMM_WORLD);
   exetime0 += MPI_Wtime();

   if(MEASURE && myproc==0)
   {
      system("/apps/power-bench/mclient -H 10.1.255.100 -e session");
      system("/apps/power-bench/mclient -H 10.1.255.100 -e log");
   }

   mapping(myproc%8, DVFS_LOW);

     
   if(myproc == SHOW)
   {
       sleep(1);
       printf("Total execution time of my_pdgemm is %.3f.\n", exetime);
       printf("Total execution time of pdgemm is %.3f.\n", exetime0);

       int i, j;
    
       printf("My PDGEMM ID AAA = %d :\n",myproc);   
        for(i=0;i<DISP_SIZE;i++)
       {
         for(j=0;j<DISP_SIZE;j++)
         	 printf("%.5lf\t", a[i*DISP_SIZE+j]);
           printf("\n");
        }
      
       printf("PDGEMM ID AAA = %d :\n",myproc);   
        for(i=0;i<DISP_SIZE;i++)
       {
         for(j=0;j<DISP_SIZE;j++)
         	 printf("%.5lf\t", A_local[i*DISP_SIZE+j]);
           printf("\n");
        }

       printf("My PDGEMM ID BBB = %d :\n",myproc);   
        for(i=0;i<DISP_SIZE;i++)
       {
         for(j=0;j<DISP_SIZE;j++)
         	 printf("%.5lf\t", b[i*DISP_SIZE+j]);
           printf("\n");
        }
      
       printf("PDGEMM ID BBB = %d :\n",myproc);   
        for(i=0;i<DISP_SIZE;i++)
       {
         for(j=0;j<DISP_SIZE;j++)
         	 printf("%.5lf\t", B_local[i*DISP_SIZE+j]);
           printf("\n");
        }

       printf("My PDGEMM ID CCC = %d :\n",myproc);   
        for(i=0;i<DISP_SIZE;i++)
       {
         for(j=0;j<DISP_SIZE;j++)
         	 printf("%.5lf\t", c[i*DISP_SIZE+j]);
           printf("\n");
        }
      
       printf("PDGEMM ID CCC = %d :\n",myproc);   
        for(i=0;i<DISP_SIZE;i++)
       {
         for(j=0;j<DISP_SIZE;j++)
         	 printf("%.5lf\t", C_local[i*DISP_SIZE+j]);
           printf("\n");
        }


   }

  
   double diff, diff_total=0.0;
   int size = MATRIX_SIZE;

   diff=diff_norm(c, C_local, size, size);
   MPI_Reduce(&diff, &diff_total, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   sleep(1);
   if(!myproc)
      printf("The total normal difference between my pdgemm  and ScaLAPACK pdgemm is %e.\n", diff_total);

   free(a); free(b); free(c);
   free(A_local);free(B_local);free(C_local);
   Cblacs_exit(1);
   /* Clean-up and close down */ 
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Comm_free(&my_row_comm); 
   MPI_Comm_free(&my_column_comm);  
   MPI_Finalize();
   return 0;
}

void Zeros(double * matrix, int r, int c){
    int i;
    for(i=0; i<(r*c) ;i++) 
       matrix[i] = 0.0;
}

void RndMatrix(double * matrix, int r, int c, int rank){
    int i;
    for(i=0; i<(r*c); i++) 
      matrix[i] = ((int)rand())%100/1000.0/(1+rank);
}

void CpyMatrix(double *a, double *b, int r, int c)
{
  int i;
  for(i=0; i<(r*c); i++)
    b[i]=a[i];
}


/*===================================================================
 *
 * Process 0 read and broadcast input data
 */
void Get_input(int p, int my_rank, int* n, int* nproc_rows, 
             int* nproc_cols,
             int* row_block_size, int* col_block_size) {

    MPI_Datatype input_datatype;

    Build_input_datatype(&input_datatype, n, nproc_rows,
            nproc_cols, row_block_size, col_block_size);
  
    
    if (my_rank == 0) 
    {
       //printf("Input 1. global size of matrix\n2. row of processors\n3. column of processors\n4. row of block size:\n5. column of block size\n");
       //scanf("%d %d %d %d %d", n, nproc_rows, 
       //       nproc_cols, row_block_size, col_block_size);
       *n=GLOBAL_M;
       *nproc_rows=PROC_NODE;
       *nproc_cols=PROC_NODE;
       *row_block_size=EACH;
       *col_block_size=EACH;

    }
  

    MPI_Bcast(n, 1, input_datatype, 0, MPI_COMM_WORLD);

    if (p < ((*nproc_rows)*(*nproc_cols))) {
        fprintf(stderr, 
           "Process %d > p = %d, nproc_rows = %d, nproc_cols = %d\n",
            my_rank, p, *nproc_rows, *nproc_cols);
        fprintf(stderr, "Process %d > Need more processes!  Quitting.", 
                my_rank);
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
}  /* Get_input */


/*===================================================================
 *
 * Build derived datatype for distributing input data
 */
void Build_input_datatype(MPI_Datatype* input_datatype, 
    int* n, int* nproc_rows, int* nproc_cols, 
    int* row_block_size, int* col_block_size) {

    int 	array_of_block_lengths[5];
    MPI_Aint 	array_of_displacements[5];
    MPI_Aint 	base_address; 
    MPI_Aint  	temp_address;
    int		i;
    MPI_Datatype  array_of_types[5]; 

    for (i = 0; i < 5; i++) {
	array_of_block_lengths[i] = 1;
        array_of_types[i] = MPI_INT;
    }
 
    /* Compute displacements from n */
    array_of_displacements[0] = 0; 
    MPI_Address(n, &base_address);
    MPI_Address(nproc_rows, &temp_address);
    array_of_displacements[1] = temp_address - base_address;
    MPI_Address(nproc_cols, &temp_address);
    array_of_displacements[2] = temp_address - base_address;
    MPI_Address(row_block_size, &temp_address);
    array_of_displacements[3] = temp_address - base_address;
    MPI_Address(col_block_size, &temp_address);
    array_of_displacements[4] = temp_address - base_address;

    MPI_Type_struct(5, array_of_block_lengths, array_of_displacements, 
                    array_of_types, input_datatype);
    MPI_Type_commit(input_datatype);

}  /* Build_input_datatype */
 

/*===================================================================
 *
 * Simplified wrapper for ScaLAPACK routine numroc                   
 * numroc computes minimum number of rows or columns needed to store 
 * a process' piece of a distributed array                          
 */
int Get_dimension(int order, int block_size, int my_process_row_or_col,
	int nproc_rows_or_cols) {

    int first_process_row_or_col = 0;  /* Assume array begins in row */
                                       /* or column 0.               */
    int return_val;

    extern int numroc_(int* order, int* block_size, 
           int* my_process_row_or_col, int* first_process_row_or_col,
           int* nproc_rows_or_cols);

    return_val = numroc_(&order, &block_size, &my_process_row_or_col,
                 &first_process_row_or_col, &nproc_rows_or_cols);
    return return_val;
}  /* Get_dimension */

/*===================================================================
 * 
 * Simplified wrapper for the ScaLAPACK routine descinit, which
 * initializes the descrip array associated to each distributed
 * matrix or vector.
 */
void Build_descrip(int my_rank, char* name, int* descrip, 
                   int m, int n, int row_block_size,
                   int col_block_size, int blacs_grid, 
                   int leading_dim) {
    int first_proc_row = 0;  /* Assume all distributed arrays begin */
    int first_proc_col = 0;  /* in process row 0, process col 0     */
    int error_info;

    extern void descinit_(int* descrip, int* m, int* n, 
                int* row_block_size, int* col_block_size, 
                int* first_proc_row, int* first_proc_col, 
                int* blacs_grid, int* leading_dim,
                int* error_info);

    descinit_(descrip, &m, &n, &row_block_size, &col_block_size,
              &first_proc_row, &first_proc_col, &blacs_grid,
              &leading_dim, &error_info);

    if (error_info != 0) {
        fprintf(stderr, "Process %d > Descinit for b failed.\n",
                my_rank);
        fprintf(stderr, "Process %d > error_info = %d\n", 
                my_rank, error_info);
        fprintf(stderr, "Process %d > Quitting.\n", my_rank);
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
}  /* Build_descrip */

void Mat_vect_mult(int m, int n, double* A_local, int* A_descrip, 
                  double* B_local, int* B_descrip, 
		  double* C_local, int* C_descrip) {
    int k = n;
    char transpose = 'N';   /* Don't take the transpose of A */
    double alpha = 1.0;
    double beta = 1.0;
    int first_row_A = 1;    /* Don't use submatrices or subvectors */
    int first_col_A = 1;    /* Multiply all of x by A              */
    int first_row_B = 1;
    int first_col_B = 1;
    int first_row_C = 1;
    int first_col_C = 1;
    
   extern void pdgemm_(char* transA, char* transB, int* m, int* n, int* k,
             double* alpha, double* A_local, int* ia, int* ja, int *A_descrip,
                            double* B_local, int* ib, int* jb, int* B_descrip,
              double* beta, double* C_local, int* ic, int* jc, int* C_descrip);
   pdgemm_(&transpose, &transpose, &m, &n, &k,
          &alpha, A_local, &first_row_A, &first_col_A, A_descrip,
                  B_local, &first_row_B, &first_col_B, B_descrip,
           &beta, C_local, &first_row_C, &first_col_C, C_descrip);

   
}  /* Mat_vect_mult */

double diff_norm(double* a, double *b, int r, int c)
{
  int i;
  double sum=0.0;
  for(i=0;i<r*c;i++)
  {
      sum += (a[i]-b[i])*(a[i]-b[i]);
  }

  return sum;
}
