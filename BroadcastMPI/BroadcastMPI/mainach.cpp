#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
int size, rank;
void random(char *a, int size)
{
   while(size--) *a++ = rand() % 100;
}
int getIteration(const int prank){
	int i;
	for(i=1;prank>=1<<i;++i);
	return i;
}
int MyBroadcastSend(void *buffer, int count, MPI_Datatype datatype, int root){
	MPI_Status status;
    int prank=(rank+(size-root))%size;
	int iter=0;
	if(root!=rank){
		iter=getIteration(prank);
		fflush(stdout);
		MPI_Recv(buffer,count,datatype, MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
	}
	int ln=1<<iter;
	while(ln+prank<size){
		MPI_Send(buffer,count,datatype,(ln+prank+root)%size,1,MPI_COMM_WORLD);
		ln<<=1;
	}
	return 0;
}
int main(int argc, char *argv[]){
	char buffer[1024*16]={0};
    int length=atoi(argv[1]);
    random(buffer,length);
	double start,end;
    const int root=0;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
    if(rank==0)start=MPI_Wtime();
    MyBroadcastSend(buffer,atoi(length),MPI_CHAR,root);
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0){
        end=MPI_Wtime();
        printf("Time my broadcast:%lf",end-start);
    }
	MPI_Finalize();
	return 0;
}
