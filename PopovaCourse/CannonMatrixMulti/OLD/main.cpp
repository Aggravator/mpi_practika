#include <mpi.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
int size, rank;
double fRand(double fMin, double fMax){
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}
void Multi(double *a,double *b,double *dest,int size){
	double *temp=new double[size];
	for(int i=0;i<size;++i){
		memset(temp,0,size*sizeof(double));
		for(int j=0;j<size;++j){
			for(int ij=0;ij<size;++ij)temp[j]+=a[i*size+ij]*b[ij*size+j];
		}
		memcpy(dest+i*size,temp,size*sizeof(double));
	}
}
void Sum(double *a,double *b,double *dest,int size){
	for(int i=0;i<size;++i)
		for(int j=0;j<size;++j)
			dest[i*size+j]=a[i*size+j]+b[i*size+j];
}
void printBlockRow(double *a,int size,int blockC){
	for(int i=0;i<size;++i){
		for(int j=0;j<blockC;++j)
			for(int ij=0;ij<size;++ij)
				printf("%lf ",a[i*size+ij+j*size*size]);
		printf("\n");
	}
}
bool isInteger(const double d,int &i){
	const static double EPS=0.00000001;
	if(fabs(ceil(d)-d)<EPS){
		i=ceil(d);
		return true;
	}
	return false;
}
int main(int argc, char *argv[]){
MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm torus;
	int matrixSize=atoi(argv[1]);
	int blockSize,cartSize;//blockSize - размер блока матрицы cartSize - размер декартовой топологии
	//Проверяем размерность топологии и блоков матрицы (обе должны быть квадратными)
	if(isInteger(sqrt((double)size),cartSize) && isInteger(sqrt(((double)matrixSize*matrixSize)/size),blockSize)){
		MPI_Request request;
		MPI_Status status;
		int zzrank;//rank процессора с адресом 0,0
		double start,end;
		{
			int dims[]={cartSize,cartSize};
			int periods[]={true,true};
			MPI_Cart_create(MPI_COMM_WORLD,2,dims,periods,true,&torus);
		}
		double *as=new double[blockSize*blockSize];//исходный блок из массива a
		double *bs=new double[blockSize*blockSize];//исходный блок из массива b
		double *t=new double[blockSize*blockSize];//переданный блок
		double *c=new double[blockSize*blockSize];//блок исходной матрицы
		memset(c,0,sizeof(double)*blockSize*blockSize);
		//Создание декартовой топологии и первичная инициализация
		int coords[2];
		{
			int nrank;
			MPI_Comm_rank(torus,&nrank);
			MPI_Cart_coords(torus,nrank,2,coords);
		}
		if(coords[0]==0 && coords[1]==0){
			int tcoords[]={0,1};
			int trank,j=1;
			for(int i=0;i<cartSize;++i){
				for(j;j<cartSize;++j){
					for(int ij=0;ij<blockSize*blockSize;++ij){
						as[ij]=fRand(-100000,100000);
						bs[ij]=fRand(-100000,100000);
					}
					tcoords[0]=i;tcoords[1]=j;
					if(i!=0)tcoords[1]-=i;
					MPI_Cart_rank(torus,tcoords,&trank);
					MPI_Send(as,blockSize*blockSize,MPI_DOUBLE,trank,1,torus);
					tcoords[1]=j;
					if(j!=0)tcoords[0]-=j;
					MPI_Cart_rank(torus,tcoords,&trank);
					MPI_Send(bs,blockSize*blockSize,MPI_DOUBLE,trank,2,torus);
				}
				j=0;
			}
			for(int ij=0;ij<blockSize*blockSize;++ij){
				as[ij]=fRand(-100000,100000);
				bs[ij]=fRand(-100000,100000);
			}
		}else{
			MPI_Irecv(as,blockSize*blockSize,MPI_DOUBLE,MPI_ANY_SOURCE,1,torus,&request);
			MPI_Recv(bs,blockSize*blockSize,MPI_DOUBLE,MPI_ANY_SOURCE,2,torus,&status);
			zzrank=status.MPI_SOURCE;
		}
		MPI_Barrier(torus);
		//Умножение
		if(coords[0]==0 && coords[1]==0) start=MPI_Wtime();
		int tcoords[]={0,0};
		int trank;
		Multi(as,bs,t,blockSize);
		Sum(c,t,c,blockSize);
		for(int i=0;i<cartSize-1;++i){
			tcoords[0]=coords[0];
			tcoords[1]=coords[1]-1;
			MPI_Cart_rank(torus,tcoords,&trank);
			MPI_Isend(as,blockSize*blockSize,MPI_DOUBLE,trank,3,torus,&request);
			MPI_Recv(t,blockSize*blockSize,MPI_DOUBLE,MPI_ANY_SOURCE,3,torus,&status);
			MPI_Wait(&request,&status);
			memcpy(as,t,blockSize*blockSize*sizeof(double));
			tcoords[0]=coords[0]-1;
			tcoords[1]=coords[1];
			MPI_Cart_rank(torus,tcoords,&trank);
			MPI_Isend(bs,blockSize*blockSize,MPI_DOUBLE,trank,4,torus,&request);
			MPI_Recv(t,blockSize*blockSize,MPI_DOUBLE,MPI_ANY_SOURCE,4,torus,&status);
			MPI_Wait(&request,&status);
			memcpy(bs,t,blockSize*blockSize*sizeof(double));
			Multi(as,bs,t,blockSize);
			Sum(c,t,c,blockSize);
		}
		//Сбор и вывод блоков результирующей матрицы
		if(coords[0]==0 && coords[1]==0){
			double *blockRow=new double[blockSize*blockSize*cartSize];
			memcpy(blockRow,c,blockSize*blockSize*sizeof(double));
			tcoords[0]=0;
			tcoords[1]=1;
			for(tcoords[0];tcoords[0]<cartSize;++tcoords[0]){
				for(tcoords[1];tcoords[1]<cartSize;++tcoords[1]){
					MPI_Cart_rank(torus,tcoords,&trank);
					MPI_Recv(blockRow+tcoords[1]*blockSize*blockSize,blockSize*blockSize,MPI_DOUBLE,trank,5,torus,&status);
				}
				//printBlockRow(blockRow,blockSize,cartSize);
				tcoords[1]=0;
			}
			end=MPI_Wtime();
			printf("\ntime:%lf",end-start);
			delete[] blockRow;
		}else{
			MPI_Send(c,blockSize*blockSize,MPI_DOUBLE,zzrank,5,torus);
		}
		delete[] as;
		delete[] bs;
		delete[] t;
		delete[] c;
	}else if(rank==0){
		printf("Size is incorrect!");
	}
	MPI_Finalize();
	return 0;
}