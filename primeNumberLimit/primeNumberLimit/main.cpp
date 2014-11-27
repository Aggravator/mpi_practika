#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
const int LIMIT=100000000;
int numCountForNode(int node,int numCount,int nodeCount){
	int res=numCount/nodeCount;
	if(node<numCount%nodeCount)++res;
	return res;
}
int getFirstNum(int node,int numCount,int nodeCount){
	int res=numCount/nodeCount;
	int temp=numCount%nodeCount;
	if(node<=temp)res=(res+1)*node;
	else res=(res+1)*temp+res*(node-temp);
	return res;
}
void getPrimeNums(bool *dest,int start,int size){
	for(int i=0;i<size;++i)dest[i]=1;
	for(int i=2;i<=(int)sqrt((double)start+size);++i){
		for(int j=0;j<size;++j){
			if(i<start+j && (start+j)%i==0)
				dest[j]=0;
		}
	}
	if(start==0){
		dest[1]=0;
		dest[0]=0;
	}
}
int size, rank;
const int partSize=3000;
int main(int argc, char *argv[]){
	MPI_Status status;
	double start,end;
	int i,j,ij,activeNodeCount;
	int n,numCount;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	if(argc<2){
		int start;
		bool *nums=new bool[partSize];
		for(int i=0;true;++i){
			start=i*partSize*size+rank*partSize;
			getPrimeNums(nums,start,partSize);
			if(rank==0){
				for(int j=0;j<partSize;++j)if(nums[j]==1)printf("%ld ",start+j);
				printf("\n");
				for(int j=1;j<size;++j){
					MPI_Recv(nums,partSize,MPI_CHAR,j,0,MPI_COMM_WORLD,&status);
					start=i*partSize*size+j*partSize;
					for(int ij=0;ij<partSize;++ij)if(nums[ij]==1)printf("%ld ",start+ij);
					printf("\n");
				}
				fflush(stdout);
			}else MPI_Send(nums,partSize,MPI_CHAR,0,0,MPI_COMM_WORLD);
			//MPI_Barrier(MPI_COMM_WORLD);
		}
	}else{
		n=atoi(argv[1]);
		activeNodeCount=n<size?n:size;
		if(rank<activeNodeCount){
			numCount=numCountForNode(rank,n,activeNodeCount);
			int firstNum=getFirstNum(rank,n,size);
			bool *nums=new bool[numCount];
			getPrimeNums(nums,firstNum,numCount);
			if(rank==0){
				for(int i=0;i<numCount;++i)if(nums[i]==1)printf("%ld ",i);
				printf("\n");
				int temp;
				int shift;
				for(int i=1;i<activeNodeCount;++i){
					temp=numCountForNode(i,n,activeNodeCount);
					MPI_Recv(nums,temp,MPI_CHAR,i,0,MPI_COMM_WORLD,&status);
					shift=getFirstNum(i,n,size);
					for(int j=0;j<temp;++j)if(nums[j]==1)printf("%ld ",shift+j);
					printf("\n");
				}
			}else MPI_Send(nums,numCount,MPI_CHAR,0,0,MPI_COMM_WORLD);
			delete[] nums;
		}
	}
	MPI_Finalize();
}