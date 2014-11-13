#include "mpi.h"
#include <stdio.h>
const int n=200;
int size, rank;
char mas[n*n+n];
int main(int argc, char *argv[]){

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Finalize();
	return 0;
}