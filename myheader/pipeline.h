#ifdef PIPELINE_H
#define PIPELINE_H


extern void MY_Bcast(double *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm, int my_rank, int pnum, int BK);

//extern void MY_Bcast2(double *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm);


#endif
