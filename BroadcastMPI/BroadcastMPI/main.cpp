#include "mpi.h"
#include <cmath>
#include <stdio.h>
const int maxBuff=200;
int size, rank;
int getPseudoRank(int root, int rank,int size){
	return (rank+(size-root))%size;
}
int getRealRank(int root,int prank,int size){
	return (prank+root)%size;
}
int getIteration(int prank){
	int lastnode=1;
	int iter=0;
	while(lastnode<prank)lastnode+=floor(pow(2.0,iter++));
	return iter;
}
int getDonorRank(int prank){
	int iter=getIteration(prank);
	return prank-int(pow(2.0,iter-1));
}
int MyBroadcastSend(void *buffer, int count, MPI_Datatype datatype, int root){
	MPI_Status status;
	int prank=getPseudoRank(root,rank,size);
	int lastnode= root==rank?1:rank;
	int iter=0;
	if(root!=rank){
		printf("I %d wait message from %d\n",rank,getRealRank(root,getDonorRank(prank),size));
		MPI_Recv(buffer,count,datatype, MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
		printf("I %d accept\n",rank);
		iter=getIteration(prank);
		lastnode+=int(pow(2.0,iter));
	}
	while(lastnode<size){
		printf("I %d send message to %d\n",rank,getRealRank(root,lastnode,size));
		MPI_Send(buffer,count,datatype,getRealRank(root,lastnode,size),MPI_ANY_TAG,MPI_COMM_WORLD);
		lastnode+=int(pow(2.0,iter++));
	}
	return 0;
}

char buffer[200]={0};
int main(int argc, char *argv[]){
	char buffer[maxBuff];
	buffer[0]='a';
	buffer[1]='d';
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	if(rank==2)	MyBroadcastSend(buffer,2,MPI_CHAR,rank);
	else MyBroadcastSend(buffer,2,MPI_CHAR,2);
	MPI_Finalize();
}