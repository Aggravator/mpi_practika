#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}
int rowCountForNode(int node,int rowCount,int nodeCount){
	int res=rowCount/nodeCount;
	if(node<rowCount%nodeCount)++res;
	return res;
}
int getIndexFirstRow(int node,int rowCount,int nodeCount){
	int res=rowCount/nodeCount;
	int temp=rowCount%nodeCount;
	if(node<=temp)res=(res+1)*node;
	else res=(res+1)*temp+res*(node-temp);
	return res;
}
int getRowOwner(int row,int rowCount,int nodeCount){
	int temp=0;
	int p=0;
	do{
		temp+=rowCountForNode(p++,rowCount,nodeCount);
	}while(row>=temp);
	return p-1;
}
const int n=7000;
int size, rank;
const double eps=0.0001;
int main(int argc, char *argv[]){
	MPI_Status status;
	int rowCount,activeNodeCount;
	int i,j,ij;
	double start,end;
	double *arrPart;
	double *mas;
	double *solution;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	activeNodeCount=size>n?n:size;
	if(rank==0){
		mas=new double[n*n+n];
		solution=new double[n];
		for(i=0;i<n;++i) solution[i]=fRand(-100,100);
		double temp;
		for(i=0;i<n;++i){
			temp=0.0;
			for(j=0;j<n;++j)
				temp+=(mas[i*(n+1)+j]=fRand(-100,100))*solution[j];
			mas[i*(n+1)+n]=temp;
		}
		rowCount=rowCountForNode(0,n,activeNodeCount);
		int shift=rowCount*n+rowCount;
		for(i=1;i<activeNodeCount;++i){
			rowCount=rowCountForNode(i,n,activeNodeCount);
			MPI_Send(mas+shift,rowCount*n+rowCount,MPI_DOUBLE,i,0,MPI_COMM_WORLD);
			shift+=rowCount*n+rowCount;
		}
		arrPart=mas;
	}
	rowCount=rowCountForNode(rank,n,activeNodeCount);
	if(rank!=0){
		arrPart=new double[n*rowCount+rowCount];
		MPI_Recv(arrPart,n*rowCount+rowCount,MPI_DOUBLE,0,0,MPI_COMM_WORLD,&status);
	}
	int rowIndex=getIndexFirstRow(rank,n,size);
	double *a=new double[n+1];
	double k;
	start=MPI_Wtime();
	//Прямой ход
	for(i=0;i<rowIndex;++i){
		MPI_Recv(a,n+1,MPI_DOUBLE,getRowOwner(i,n,activeNodeCount),MPI_ANY_TAG,MPI_COMM_WORLD,&status);
		for(j=0;j<rowCount;++j){
			k=-arrPart[j*(n+1)+i]/a[i];
			for(ij=i;ij<=n;++ij)arrPart[j*(n+1)+ij]+=a[ij]*k;
		}
	}
	for(i=0;i<rowCount;++i){
		for(j=rank+1;j<activeNodeCount;++j)MPI_Send(arrPart+i*(n+1),n+1,MPI_DOUBLE,j,0,MPI_COMM_WORLD);
		for(j=i+1;j<rowCount;++j){
			k=-arrPart[j*(n+1)+rowIndex+i]/arrPart[i*(n+1)+rowIndex+i];
			for(ij=rowIndex+i;ij<=n;++ij)arrPart[j*(n+1)+ij]+=arrPart[i*(n+1)+ij]*k;
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	//Обратный ход
	for(i=n-1;i>=rowIndex+rowCount;--i){
		MPI_Recv(a,n+1,MPI_DOUBLE,getRowOwner(i,n,activeNodeCount),MPI_ANY_TAG,MPI_COMM_WORLD,&status);
		for(j=rowCount-1;j>=0;--j){
			k=-arrPart[j*(n+1)+i]/a[i];
			arrPart[j*(n+1)+i]+=a[i]*k;
			arrPart[j*(n+1)+n]+=a[n]*k;
		}
	}
	for(i=rowCount-1;i>=0;--i){
		for(j=rank-1;j>=0;--j)MPI_Send(arrPart+i*(n+1),n+1,MPI_DOUBLE,j,0,MPI_COMM_WORLD);
		for(j=i-1;j>=0;--j){
			k=-arrPart[j*(n+1)+rowIndex+i]/arrPart[i*(n+1)+rowIndex+i];
			arrPart[j*(n+1)+rowIndex+i]+=arrPart[i*(n+1)+rowIndex+i]*k;
			arrPart[j*(n+1)+n]+=arrPart[i*(n+1)+n]*k;
		}
	}
	end=MPI_Wtime();
	bool result=true;
	if(rank==0){
		for(i=0;i<rowCount;++i)
			if(fabs(arrPart[(n+1)*(i+1)-1]/arrPart[(n+1)*i+i]-solution[i])>eps)result=false;
		for(j=1;j<activeNodeCount;++j){
			MPI_Recv(a,rowCountForNode(j,n,activeNodeCount),MPI_DOUBLE,j,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
			for(ij=0;ij<rowCountForNode(j,n,activeNodeCount);++ij)if(fabs(a[ij]-solution[i++])>eps)result=false;
		}
		printf("%s\n",result?"Сошлось":"Печалька");
		printf("Время:%lf",end-start);
		delete[] mas;
		delete[] solution; 
	}else{
		for(i=0;i<rowCount;++i)arrPart[i]=arrPart[(i+1)*(n+1)-1]/arrPart[rowIndex+i*(n+1)+i];
		MPI_Send(arrPart,rowCount,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
		delete[] arrPart;
	}
	MPI_Finalize();
	return 0;
}