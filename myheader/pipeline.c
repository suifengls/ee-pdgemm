#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include <time.h>
#include "malloc.h"

#include "pipeline.h"

void MY_Bcast(double *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm, int my_rank, int pnum, int BK)
{
  //  double time;
   int i;
   int numbk;
   int last=0;
   //   int count = row*column;
   //int first = root%2+1;
   MPI_Status status;

   numbk=count/BK;
   last=count%BK;
   //printf("num_bk=%d, last_size = %d\n",numbk, last);

   //MPI_Barrier(MPI_COMM_WORLD);
   //time=-MPI_Wtime();

   if(my_rank==root)
   {
      if(last==0)
      {
         for(i=0;i<numbk;i++)
            MPI_Send(&buffer[i*BK], BK, datatype, (root+1)%pnum, 421, comm);
      }
      else
      {
         for(i=0;i<numbk;i++)
            MPI_Send(&buffer[i*BK], BK, datatype, (root+1)%pnum, 421, comm);

         MPI_Send(&buffer[(numbk)*BK], last, datatype, (root+1)%pnum, 628, comm);
      }
   }
   else if(my_rank==(root-1+pnum)%pnum)
   {
      if(last==0)
      {
         for(i=0;i<numbk;i++)
            MPI_Recv(&buffer[i*BK], BK, datatype, (root-2+pnum)%pnum, 421, comm, &status);
      }
      else
      {
         for(i=0;i<numbk;i++)
            MPI_Recv(&buffer[i*BK], BK, datatype, (root-2+pnum)%pnum, 421, comm, &status);

         MPI_Recv(&buffer[(numbk)*BK], last, datatype, (root-2+pnum)%pnum, 628, comm, &status);
      }
   }
   else if( (my_rank != root)&&(my_rank != (root-1+pnum)%pnum) )
   {
      if(last==0)
      {
         for(i=0;i<numbk;i++)
         {
            MPI_Recv(&buffer[i*BK], BK, datatype, (my_rank-1+pnum)%pnum, 421, comm, &status);
            MPI_Send(&buffer[i*BK], BK, datatype, (my_rank+1)%pnum, 421, comm);
         }
       }
       else
       {
         for(i=0;i<numbk;i++)
         {
            MPI_Recv(&buffer[i*BK], BK, datatype, (my_rank-1+pnum)%pnum, 421, comm, &status);
            MPI_Send(&buffer[i*BK], BK, datatype, (my_rank+1)%pnum, 421, comm);
         }
            MPI_Recv(&buffer[(numbk)*BK], last, datatype, (my_rank-1+pnum)%pnum, 628, comm, &status);
            MPI_Send(&buffer[(numbk)*BK], last, datatype, (my_rank+1)%pnum, 628, comm);
       }
   }

   //MPI_Barrier(MPI_COMM_WORLD);
   //time+=MPI_Wtime();
   //if(my_rank==0)
   //    printf("The execution of MY_Bcast() is : %.6lf\n",time);
}

/*
void MY_Bcast2(double *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm)
{
   //double time
   int BK = 1024;
   int i;
   int numbk;
   int last=0;
   MPI_Status status;
   int my_rank;
   int pnum;
 
   MPI_Comm_size(comm, &pnum);
   MPI_Comm_rank(comm, &my_rank);

   numbk=count/BK;
   last=count%BK;
   //   printf("rank=%d, pnum =%d, root = %d\n", my_rank, pnum, root);

   //MPI_Barrier(MPI_COMM_WORLD);
   //time=-MPI_Wtime();

   if(my_rank==root)
   {
      if(last==0)
      {
         for(i=0;i<numbk;i++)
            MPI_Send(&buffer[i*BK], BK, datatype, (root+1)%pnum, 421, comm);
      }
      else
      {
         for(i=0;i<numbk-1;i++)
            MPI_Send(&buffer[i*BK], BK, datatype, (root+1)%pnum, 421, comm);

         MPI_Send(&buffer[(numbk)*BK], last, datatype, (root+1)%pnum, 421, comm);
      }
   }
   else if(my_rank==(root-1+pnum)%pnum)
   {
      if(last==0)
      {
         for(i=0;i<numbk;i++)
            MPI_Recv(&buffer[i*BK], BK, datatype, (root-2+pnum)%pnum, 421, comm, &status);
      }
      else
      {
         for(i=0;i<numbk-1;i++)
            MPI_Recv(&buffer[i*BK], BK, datatype, (root-2+pnum)%pnum, 421, comm, &status);

         MPI_Recv(&buffer[(numbk)*BK], last, datatype, (root-2+pnum)%pnum, 421, comm, &status);
      }
   }
   else //if(my_rank%2==first)
   {
      if(last==0)
      {
         for(i=0;i<numbk;i++)
         {
            MPI_Recv(&buffer[i*BK], BK, datatype, (my_rank-1+pnum)%pnum, 421, comm, &status);
            MPI_Send(&buffer[i*BK], BK, datatype, (my_rank+1)%pnum, 421, comm);
         }
       }
       else
       {
         for(i=0;i<numbk-1;i++)
         {
            MPI_Recv(&buffer[i*BK], BK, datatype, (my_rank-1+pnum)%pnum, 421, comm, &status);
            MPI_Send(&buffer[i*BK], BK, datatype, (my_rank+1)%pnum, 421, comm);
         }
            MPI_Recv(&buffer[(numbk)*BK], last, datatype, (my_rank-1+pnum)%pnum, 421, comm, &status);
            MPI_Send(&buffer[(numbk)*BK], last, datatype, (my_rank+1)%pnum, 421, comm);
       }
   }

   //MPI_Barrier(MPI_COMM_WORLD);
   //time+=MPI_Wtime();
   //if(my_rank==0)
   //    printf("The execution of MY_Bcast() is : %.6lf\n",time);
}
*/
