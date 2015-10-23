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
	double start,end;
	const int root=0;
	int dataSize[]={16,256,512,1024,1024*4,1024*8,1024*16};
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	for(int i=0;i<sizeof(dataSize)/4;++i){
		random(buffer,dataSize[i]);
		if(rank==0)start=MPI_Wtime();
		MyBroadcastSend(buffer,dataSize[i],MPI_CHAR,root);
		MPI_Barrier(MPI_COMM_WORLD);
		if(rank==0){
			end=MPI_Wtime();
			printf("\nData size:%d\n",dataSize[i]);
			printf("	Time my broadcast:%lf",end-start);
			random(buffer,dataSize[i]);
			start=MPI_Wtime();
		}
		MPI_Bcast(buffer,dataSize[i],MPI_CHAR,root,MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		if(rank==0){
			end=MPI_Wtime();
			printf("\n	Time standart broadcast:%lf",end-start);
		}
	}	
	MPI_Finalize();
	return 0;
}
