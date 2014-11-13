#include "mpi.h"
#include <stdio.h>
const int maxBuff=200;
int size, rank;
int getIteration(const int prank){
	int lastnode=1,iter=1;
	for(lastnode;lastnode<prank;lastnode<<=1)++iter;
	return iter;
}
int MyBroadcastSend(void *buffer, int count, MPI_Datatype datatype, int root){
	MPI_Status status;
	int prank=(rank+(size-root))%size;
	int iter=0;
	if(root!=rank){
		iter=getIteration(prank);
		MPI_Recv(buffer,count,datatype, MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
	}
	int lastnode=1<<iter;
	while(lastnode+prank<size){
		MPI_Send(buffer,count,datatype,(lastnode+prank+root)%size,1,MPI_COMM_WORLD);
		lastnode<<=1;
	}
	return 0;
}
int main(int argc, char *argv[]){
	char buffer[maxBuff];
	const int root=4;
	buffer[0]='a';
	buffer[1]='d';
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MyBroadcastSend(buffer,2,MPI_CHAR,root);
	MPI_Finalize();
	return 0;
}