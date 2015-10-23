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
	double totalTime=MPI_Wtime(),readingTime,writingTime,perfomanceTime;
	int matrixSize=atoi(argv[1]);
	int blockSize,cartSize;//blockSize - ðàçìåð áëîêà ìàòðèöû cartSize - ðàçìåð äåêàðòîâîé òîïîëîãèè
	//Ïðîâåðÿåì ðàçìåðíîñòü òîïîëîãèè è áëîêîâ ìàòðèöû (îáå äîëæíû áûòü êâàäðàòíûìè)

	MPI_File f;
	MPI_Status s;
	if(isInteger(sqrt((double)size),cartSize) && isInteger(sqrt(((double)matrixSize*matrixSize)/size),blockSize)){
		MPI_Request request;
		MPI_Status status;
		int zzrank=0;//rank ïðîöåññîðà ñ àäðåñîì 0,0
		double start,end;
		int coords[2];
		{
			int dims[]={cartSize,cartSize};
			int periods[]={true,true};
			MPI_Cart_create(MPI_COMM_WORLD,2,dims,periods,true,&torus);
			int nrank;
			MPI_Comm_rank(torus,&nrank);
			MPI_Cart_coords(torus,nrank,2,coords);
		}
		double *as=new double[blockSize*blockSize];//èñõîäíûé áëîê èç ìàññèâà a
		double *bs=new double[blockSize*blockSize];//èñõîäíûé áëîê èç ìàññèâà b
		double *t=new double[blockSize*blockSize];//ïåðåäàííûé áëîê
		double *c=new double[blockSize*blockSize];//áëîê èñõîäíîé ìàòðèöû
		memset(c,0,sizeof(double)*blockSize*blockSize);
		//Ñîçäàíèå äåêàðòîâîé òîïîëîãèè è ïåðâè÷íàÿ èíèöèàëèçàöèÿ
		MPI_Datatype cOutType;
		//Input matrixes A and B
		{
			MPI_Datatype mysubarray;
			int subArray[] = {blockSize,blockSize};
			int bigsizes[2] = { matrixSize, matrixSize };
			int starts[] = {coords[0]*blockSize,coords[1]*blockSize};
			MPI_Type_create_subarray(2,bigsizes,subArray,starts,MPI_ORDER_C, MPI_DOUBLE,&cOutType);
			starts[1]=((coords[1]+coords[0])%cartSize)*blockSize;
			MPI_Type_create_subarray(2, bigsizes, subArray, starts, MPI_ORDER_C, MPI_DOUBLE, &mysubarray);
			readingTime=MPI_Wtime();
			MPI_File_open(MPI_COMM_WORLD, "A", MPI_MODE_RDONLY, MPI_INFO_NULL, &f);
			MPI_File_set_view(f, 0, MPI_DOUBLE, mysubarray, "native", MPI_INFO_NULL);
			if (MPI_File_read(f, as, blockSize*blockSize, MPI_DOUBLE, &s) != MPI_SUCCESS) printf("OOOOOPS A READING!\n");
			MPI_File_close(&f);
			MPI_Type_free(&mysubarray);
			starts[0]=((coords[0]+coords[1])%cartSize)*blockSize;
			starts[1]=coords[1]*blockSize;
			MPI_Type_create_subarray(2, bigsizes, subArray, starts, MPI_ORDER_C, MPI_DOUBLE, &mysubarray);
			MPI_File_open(MPI_COMM_WORLD, "B", MPI_MODE_RDONLY, MPI_INFO_NULL, &f);
			MPI_File_set_view(f, 0, MPI_DOUBLE, mysubarray, "native", MPI_INFO_NULL);
			if (MPI_File_read(f, bs, blockSize*blockSize, MPI_DOUBLE, &s) != MPI_SUCCESS) printf("OOOOOPS B READING!\n");
			MPI_File_close(&f);
			MPI_Type_free(&mysubarray);
			readingTime=MPI_Wtime()-readingTime;
		}

		MPI_Barrier(torus);
		//Óìíîæåíèå
		perfomanceTime=MPI_Wtime();
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
		perfomanceTime=MPI_Wtime()-perfomanceTime;
		//Output evaluated matrix
		{
			writingTime=MPI_Wtime();
			MPI_File_open(MPI_COMM_WORLD, "C", MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &f);
			MPI_File_set_view(f, 0, MPI_DOUBLE, cOutType, "native", MPI_INFO_NULL);
			MPI_File_write(f, c, blockSize*blockSize, MPI_DOUBLE, &s);
			MPI_File_close(&f);
			writingTime=MPI_Wtime()-writingTime;
		}
		if(rank==0){
			printf("Execution time:%f\n",MPI_Wtime()-totalTime);
			printf("Reading time:%f\n",readingTime);
			printf("Performance time:%f\n",perfomanceTime);
			printf("Writing time:%f\n",writingTime);
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

