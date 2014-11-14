#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
float fRand(float fMin, float fMax){
    float f = (float)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}
int rowCountForNode(int node,int rowCount,int nodeCount){
	int res=rowCount/nodeCount;
	if(node<rowCount%nodeCount)++res;
	return res;
}
int size, rank;
int main(int argc, char *argv[]){
	MPI_Status status;
	float *matrix,*vector,*answer;
	float start,end;
	int n,i,j,ij,activeNodeCount,rowCount;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	n=atoi(argv[1]);
	activeNodeCount=n<size?n:size;
	rowCount=rowCountForNode(rank,n,activeNodeCount);
	if(rank<activeNodeCount){
		matrix=new float[rowCount*n];
		vector=new float[n];
		::fflush(stdout);
		if(rank==0){
			int temp;
			for(i=0;i<n;++i)vector[i]=fRand(-100.0,100.0);
			for(i=1;i<activeNodeCount;++i){
				temp=rowCountForNode(i,n,activeNodeCount);
				for(j=0;j<temp;++j){
					for(ij=0;ij<n;++ij){
						matrix[j*n+ij]=fRand(-100.0,100.0);
					}
				}
				MPI_Send(vector,n,MPI_FLOAT,i,0,MPI_COMM_WORLD);
				MPI_Send(matrix,temp*n,MPI_FLOAT,i,1,MPI_COMM_WORLD);
			}
			for(i=0;i<rowCount;++i)for(j=0;j<n;++j)matrix[i*n+j]=fRand(-100.0,100.0);
			answer=new float[n];
		}else{
			MPI_Recv(vector,n,MPI_FLOAT,0,0,MPI_COMM_WORLD,&status);
			MPI_Recv(matrix,rowCount*n,MPI_FLOAT,0,1,MPI_COMM_WORLD,&status);
			answer=new float[rowCount];
		}
		MPI_Barrier(MPI_COMM_WORLD);
		if(rank==0)start=MPI_Wtime();
		for(i=0;i<rowCount;++i)
			for(j=0;j<n;++j)answer[i]+=matrix[i*n+j]*vector[j];
		if(rank==0){
			int shift=rowCount;
			int temp;
			for(i=1;i<activeNodeCount;++i){
				temp=rowCountForNode(i,n,activeNodeCount);
				MPI_Recv(answer+shift,temp,MPI_FLOAT,i,0,MPI_COMM_WORLD,&status);
				shift+=temp;
			}
			end=MPI_Wtime();
			printf("time:%lf\n\n",end-start);
			for(i=0;i<n;++i)printf("%lf\n",answer[i]);
		}else MPI_Send(answer,rowCount,MPI_FLOAT,0,0,MPI_COMM_WORLD);
	}
	MPI_Finalize();
	return 0;
}