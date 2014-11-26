#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
double fRand(double fMin, double fMax){
    double f = (double)rand() / RAND_MAX;
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
	double *matrix,*vector,*answer;
	double start,end;
	int n,i,j,ij,activeNodeCount,rowCount;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	n=atoi(argv[1]);
	activeNodeCount=n<size?n:size;
	rowCount=rowCountForNode(rank,n,activeNodeCount);
	if(rank<activeNodeCount){
		matrix=new double[rowCount*n];
		vector=new double[n];
		if(rank==0){
			int tempRow;
			for(i=0;i<n;++i)vector[i]=fRand(-100.0,100.0);
			for(i=1;i<activeNodeCount;++i){
				tempRow=rowCountForNode(i,n,activeNodeCount);
				for(j=0;j<tempRow;++j){
					for(ij=0;ij<n;++ij){
						matrix[j*n+ij]=fRand(-100.0,100.0);
					}
				}
				MPI_Send(vector,n,MPI_DOUBLE,i,0,MPI_COMM_WORLD);
				MPI_Send(matrix,tempRow*n,MPI_DOUBLE,i,1,MPI_COMM_WORLD);
			}
			for(i=0;i<rowCount;++i)for(j=0;j<n;++j)matrix[i*n+j]=fRand(-100.0,100.0);
			answer=new double[n];
		}else{
			MPI_Recv(vector,n,MPI_DOUBLE,0,0,MPI_COMM_WORLD,&status);
			MPI_Recv(matrix,rowCount*n,MPI_DOUBLE,0,1,MPI_COMM_WORLD,&status);
			answer=new double[rowCount];
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
				MPI_Recv(answer+shift,temp,MPI_DOUBLE,i,0,MPI_COMM_WORLD,&status);
				shift+=temp;
			}
			end=MPI_Wtime();
			printf("time:%lf\n\n",end-start);
			for(i=0;i<n;++i)printf("%lf\n",answer[i]);
		}else MPI_Send(answer,rowCount,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
	}
	MPI_Finalize();
	return 0;
}