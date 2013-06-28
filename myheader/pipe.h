#ifdef PIPE_H
#define PIPE_H

extern void MY_Bcast3(double *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm);

#endif
