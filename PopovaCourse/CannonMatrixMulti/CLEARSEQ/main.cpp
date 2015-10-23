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
	int blockSize,cartSize;//blockSize - ðàçìåð áëîêà ìàòðèöû cartSize - ðàçìåð äåêàðòîâîé òîïîëîãèè
	//Ïðîâåðÿåì ðàçìåðíîñòü òîïîëîãèè è áëîêîâ ìàòðèöû (îáå äîëæíû áûòü êâàäðàòíûìè)
	double totalTime=MPI_Wtime(),readingTime,writingTime,perfomanceTime;

	if(isInteger(sqrt((double)size),cartSize) && isInteger(sqrt(((double)matrixSize*matrixSize)/size),blockSize)){
		MPI_Request request;
		MPI_Status status;
		int zzrank;//rank ïðîöåññîðà ñ àäðåñîì 0,0
		double start,end;
		{
			int dims[]={cartSize,cartSize};
			int periods[]={true,true};
			MPI_Cart_create(MPI_COMM_WORLD,2,dims,periods,true,&torus);
		}
		double *as=new double[blockSize*blockSize];//èñõîäíûé áëîê èç ìàññèâà a
		double *bs=new double[blockSize*blockSize];//èñõîäíûé áëîê èç ìàññèâà b
		double *t=new double[blockSize*blockSize];//ïåðåäàííûé áëîê
		double *c=new double[blockSize*blockSize];//áëîê èñõîäíîé ìàòðèöû
		memset(c,0,sizeof(double)*blockSize*blockSize);
		//Ñîçäàíèå äåêàðòîâîé òîïîëîãèè è ïåðâè÷íàÿ èíèöèàëèçàöèÿ
		int coords[2];
		int tcoords[]={0,0};
		int trank;
		{
			int nrank;
			MPI_Comm_rank(torus,&nrank);
			MPI_Cart_coords(torus,nrank,2,coords);
		}
		//Input matrixes A and B
		{
			readingTime=MPI_Wtime();
			if(coords[0]==0 || coords[1]==0){
				int offset=coords[0]*blockSize*blockSize*cartSize+((coords[1]+coords[0])%cartSize)*blockSize;
				double *row=as;
				FILE *ff=fopen("A","rb");
				for(int i=0;i<blockSize;++i){
					fseek(ff,(offset+i*matrixSize)*sizeof(as[0]),SEEK_SET);
					row+=fread(row,sizeof(as[0]),blockSize,ff);
				}
				row=t;
				tcoords[1]=1;
				for(tcoords[0]=0;tcoords[0]<cartSize;++tcoords[0]){
					for(tcoords[1];tcoords[1]<cartSize;++tcoords[1]){
						offset=tcoords[0]*blockSize*blockSize*cartSize+((tcoords[1]+tcoords[0])%cartSize)*blockSize;
						for(int ij=0;ij<blockSize;++ij){
							fseek(ff,(offset+ij*matrixSize)*sizeof(as[0]),SEEK_SET);
							row+=fread(row,sizeof(as[0]),blockSize,ff);
						}
						MPI_Cart_rank(torus,tcoords,&trank);
						MPI_Send(t,blockSize*blockSize,MPI_DOUBLE,trank,2,torus);
					}
					tcoords[1]=0;
				}
				fclose(ff);
				
				row=bs;
				ff=fopen("B","rb");
				offset=((coords[1]+coords[0])%cartSize)*blockSize*blockSize*cartSize+coords[1]*blockSize;
				for(int i=0;i<blockSize;++i){
					fseek(ff,(offset+i*matrixSize)*sizeof(as[0]),SEEK_SET);
					row+=fread(row,sizeof(as[0]),blockSize,ff);
				}
				row=t;
				tcoords[1]=1;
				for(tcoords[0]=0;tcoords[0]<cartSize;++tcoords[0]){
					for(tcoords[1];tcoords[1]<cartSize;++tcoords[1]){
						offset=((tcoords[1]+tcoords[0])%cartSize)*blockSize*blockSize*cartSize+tcoords[1]*blockSize;
						for(int ij=0;ij<blockSize;++ij){
							fseek(ff,(offset+ij*matrixSize)*sizeof(as[0]),SEEK_SET);
							row+=fread(row,sizeof(as[0]),blockSize,ff);
						}
						MPI_Cart_rank(torus,tcoords,&trank);
						MPI_Send(t,blockSize*blockSize,MPI_DOUBLE,trank,2,torus);
					}
					tcoords[1]=0;
				}
				fclose(ff);
			}else{
				MPI_Recv(as,blockSize*blockSize,MPI_DOUBLE,MPI_ANY_SOURCE,2,torus,&status);
				MPI_Recv(bs,blockSize*blockSize,MPI_DOUBLE,MPI_ANY_SOURCE,2,torus,&status);
			}
			readingTime=MPI_Wtime()-readingTime;
		}
		MPI_Barrier(torus);

		perfomanceTime=MPI_Wtime();
		if(coords[0]==0 && coords[1]==0) start=MPI_Wtime();
		
		
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
			if(coords[0]==0 && coords[1]==0){
				int offset=0;
				FILE *ff=fopen("C","w");
				fseek(ff,sizeof(c[0])*matrixSize*matrixSize-1,SEEK_SET);
				fwrite("\0",1,1,ff);
				double *pos=c;
				for(int i=0;i<blockSize;++i){
					fseek(ff,sizeof(c[0])*(offset+i*matrixSize),SEEK_SET);
					pos+=fwrite(pos,sizeof(c[0]),blockSize,ff);
				}
				tcoords[0]=0;
				tcoords[1]=1;
				for(tcoords[0];tcoords[0]<cartSize;++tcoords[0]){
					for(tcoords[1];tcoords[1]<cartSize;++tcoords[1]){
						MPI_Cart_rank(torus,tcoords,&trank);
						MPI_Recv(c,blockSize*blockSize,MPI_DOUBLE,trank,5,torus,&status);
						offset=(tcoords[0]*cartSize)*(blockSize*blockSize)+tcoords[1]*blockSize;
						pos=c;
						for(int i=0;i<blockSize;++i){
							fseek(ff,sizeof(c[0])*(offset+i*matrixSize),SEEK_SET);
							pos+=fwrite(pos,sizeof(c[0]),blockSize,ff);
						}
					}
					tcoords[1]=0;
				}
				fclose(ff);
			}else{
				MPI_Send(c,blockSize*blockSize,MPI_DOUBLE,zzrank,5,torus);
			}
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

