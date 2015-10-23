#include <mpi.h>
#include <omp.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
const double EPS=0.00000001;
int size, rank;
double fRand(double fMin, double fMax){
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}
void Multi(double *a,double *b,double *dest,int size){
	#pragma omp parallel shared(a,b,dest)
	{
		double *temp=new double[size];
		#pragma omp for nowait
		for(int i=0;i<size;++i){
			memset(temp,0,size*sizeof(double));
			for(int j=0;j<size;++j){
				for(int ij=0;ij<size;++ij)temp[j]+=a[i*size+ij]*b[ij*size+j];
			}
			memcpy(dest+i*size,temp,size*sizeof(double));
		}
		delete[] temp;
	}
}
void Sum(double *a,double *b,double *dest,int size){
	#pragma omp parallel shared(a,b,dest)
	{
		#pragma omp for
		for(int i=0;i<size;++i)
			for(int j=0;j<size;++j)
				dest[i*size+j]=a[i*size+j]+b[i*size+j];
	}
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
	double totalTime=MPI_Wtime(),readingTime,writingTime,perfomanceTime;
	MPI_Comm torus;
	int matrixSize=atoi(argv[1]);
	int blockSize,cartSize;//blockSize - ðàçìåð áëîêà ìàòðèöû cartSize - ðàçìåð äåêàðòîâîé òîïîëîãèè
	//Ïðîâåðÿåì ðàçìåðíîñòü òîïîëîãèè è áëîêîâ ìàòðèöû (îáå äîëæíû áûòü êâàäðàòíûìè)

	MPI_Status s;
	if(isInteger(sqrt((double)size),cartSize) && isInteger(sqrt(((double)matrixSize*matrixSize)/size),blockSize)){
		if(rank==0)printf("Cart size and matrix size is OK!\n");
		omp_set_num_threads(4);
		MPI_Request request;
		MPI_Status status;
		int zzrank;//rank ïðîöåññîðà ñ àäðåñîì 0,0
		double start,end;
		{
			int dims[]={cartSize,cartSize};
			int periods[]={true,false};
			MPI_Cart_create(MPI_COMM_WORLD,2,dims,periods,true,&torus);
		}
		double *as=new double[blockSize*blockSize];//èñõîäíûé áëîê èç ìàññèâà a
		double *bs=new double[blockSize*blockSize];//èñõîäíûé áëîê èç ìàññèâà b
		double *at=new double[blockSize*blockSize];//ïåðåäàííûé áëîê èç ìàññèâà a
		double *bt=new double[blockSize*blockSize];//ïåðåäàííûé áëîê èç ìàññèâà b
		double *c=new double[blockSize*blockSize];//áëîê èñõîäíîé ìàòðèöû
		memset(c,0,sizeof(double)*blockSize*blockSize);
		//Ñîçäàíèå äåêàðòîâîé òîïîëîãèè è ïåðâè÷íàÿ èíèöèàëèçàöèÿ
		int coords[2];
		{
			int nrank;
			MPI_Comm_rank(torus,&nrank);
			MPI_Cart_coords(torus,nrank,2,coords);
		}
		//Input matrixes A and B
		{
			readingTime=MPI_Wtime();
			int offset=(rank/cartSize)*(blockSize*blockSize*cartSize)+(rank%cartSize)*blockSize;
			double *row=as;
			FILE *ff=fopen("A","rb");
			for(int i=0;i<blockSize;++i){
				fseek(ff,(offset+i*matrixSize)*sizeof(as[0]),SEEK_SET);
				row+=fread(row,sizeof(as[0]),blockSize,ff);
			}
			fclose(ff);
			row=bs;
			ff=fopen("B","rb");
			for(int i=0;i<blockSize;++i){
				fseek(ff,(offset+i*matrixSize)*sizeof(as[0]),SEEK_SET);
				row+=fread(row,sizeof(as[0]),blockSize,ff);
			}
			fclose(ff);
			readingTime=MPI_Wtime()-readingTime;
		}
		if(rank==0)printf("A and B array have been read!\n");
		MPI_Barrier(torus);
		//Óìíîæåíèå
		if(coords[0]==0 && coords[1]==0) start=MPI_Wtime();
		int tcoords[]={0,0};
		int trank;
		perfomanceTime=MPI_Wtime();
		for(int i=0;i<cartSize;++i){
			if(coords[1]==(coords[0]+i)%cartSize){
				tcoords[0]=coords[0];
				for(int j=0;j<cartSize;++j){
					if(j==coords[1])continue;
					tcoords[1]=j;
					MPI_Cart_rank(torus,tcoords,&trank);
					MPI_Isend(as,blockSize*blockSize,MPI_DOUBLE,trank,3,torus,&request);
				}
				memcpy(at,as,blockSize*blockSize*sizeof(double));
			}else{
				MPI_Recv(at,blockSize*blockSize,MPI_DOUBLE,MPI_ANY_SOURCE,3,torus,&status);
			}
			if(i==0){
				memcpy(bt,bs,blockSize*blockSize*sizeof(double));
			}else{
				tcoords[0]=coords[0]-1;
				tcoords[1]=coords[1];
				MPI_Cart_rank(torus,tcoords,&trank);
				MPI_Isend(bs,blockSize*blockSize,MPI_DOUBLE,trank,4,torus,&request);
				MPI_Recv(bt,blockSize*blockSize,MPI_DOUBLE,MPI_ANY_SOURCE,4,torus,&status);
			}
			Multi(at,bt,at,blockSize);
			Sum(c,at,c,blockSize);
		}
		perfomanceTime=MPI_Wtime()-perfomanceTime;
		if(rank==0)printf("Multiplication has been performed\n");
		//Output evaluated matrix
		{
			writingTime=MPI_Wtime();
			if(rank==0){
				FILE *ff=fopen("C","w");
				fseek(ff,sizeof(c[0])*matrixSize*matrixSize-1,SEEK_SET);
				printf("C File size %d\n",sizeof(c[0])*matrixSize*matrixSize);
				fwrite("\0",1,1,ff);
				fclose(ff);
			}
			MPI_Barrier(MPI_COMM_WORLD);
			int offset=(rank/cartSize)*(blockSize*blockSize*cartSize)+(rank%cartSize)*blockSize;
			double *row=c;
			FILE *ff=fopen("C","r+");
			fseek(ff,offset,SEEK_SET);
			for(int i=0;i<blockSize;++i){
				fseek(ff,(offset+i*matrixSize)*sizeof(c[0]),SEEK_SET);
				fwrite(c,sizeof(c[0]),blockSize,ff);
			}
			fclose(ff);
			writingTime=MPI_Wtime()-writingTime;
		}
		if(rank==0)printf("C array has been written\n");
		if(rank==0){
			printf("Execution time:%f\n",MPI_Wtime()-totalTime);
			printf("Reading time:%f\n",readingTime);
			printf("Performance time:%f\n",perfomanceTime);
			printf("Writing time:%f\n",writingTime);
		}
		delete[] as;
		delete[] bs;
		delete[] at;
		delete[] bt;
		delete[] c;
	}else if(rank==0){
		printf("Size is incorrect!");
	}
	MPI_Finalize();
	return 0;
}

