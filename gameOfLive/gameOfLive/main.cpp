#include <mpi.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
double fRand(double fMin, double fMax){
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}
bool isInteger(const double d,int &i){
	const static double EPS=0.00000001;
	if(fabs(ceil(d)-d)<EPS){
		i=ceil(d);
		return true;
	}
	return false;
}
void generateBlock(char *block,int size){
	for(int i=0;i<size;++i)
		for(int j=0;j<size;++j)
			block[i*size+j]=fRand(0,3)>2 ? 1 :0;
}
void printBlock(char *block,int size){
	for(int i=0;i<size;++i){
		for(int j=0;j<size;++j){
			printf("%d ",block[i*size+j]);
		}
		printf("\n");
	}
}
void blockToLay(char* block,char** lay,char *lc,char *rc,int size){
	for(int i=0;i<size;++i){
		for(int j=0;j<size;++j)
			lay[1+i][1+j]=block[i*size+j];
		lc[i]=lay[1+i][1];
		rc[i]=lay[1+i][size];
	}
}
void ghostShare(MPI_Comm torus,char** lay,char *lc,char *rc,int coords[2],int bS,int cS){
	int tcoords[]={coords[0],coords[1]};
	int trank;
	static char *temp=new char[bS];
	MPI_Status status;
	MPI_Request request;
	/*tcoords[1]=coords[1]-1;
	MPI_Cart_rank(torus,tcoords,&trank);
	if(coords[1]%2==0)MPI_Sendrecv(lc,bS,MPI_CHAR,trank,5,temp,bS,MPI_CHAR,trank,6,torus,&status);
	else MPI_Sendrecv(rc,bS,MPI_CHAR,trank,6,temp,bS,MPI_CHAR,trank,5,torus,&status);
	for(int i=0;i<bS;++i)lay[i][bS+1]=temp[i];*/
	tcoords[1]=coords[1]-1;
	MPI_Cart_rank(torus,tcoords,&trank);
	MPI_Isend(lc,bS,MPI_CHAR,trank,5,torus,&request);
	MPI_Recv(temp,bS,MPI_CHAR,MPI_ANY_SOURCE,5,torus,&status);
	for(int i=0;i<bS;++i)lay[i+1][bS+1]=temp[i];
	MPI_Wait(&request,&status);

	tcoords[1]=coords[1]+1;
	MPI_Cart_rank(torus,tcoords,&trank);
	MPI_Isend(rc,bS,MPI_CHAR,trank,6,torus,&request);
	MPI_Recv(temp,bS,MPI_CHAR,MPI_ANY_SOURCE,6,torus,&status);
	for(int i=0;i<bS;++i)lay[i+1][0]=temp[i];
	MPI_Wait(&request,&status);

	tcoords[1]=coords[1];
	tcoords[0]=coords[0]-1;
	printf("%d\n",lay[0][0]);
	MPI_Cart_rank(torus,tcoords,&trank);
	MPI_Isend(lay[0],bS+2,MPI_CHAR,trank,7,torus,&request);
	printf("bot gostsent!");
	fflush(stdout);
	MPI_Recv(lay[bS+2],bS+2,MPI_CHAR,MPI_ANY_SOURCE,7,torus,&status);
	MPI_Wait(&request,&status);
	printf("bot gostsent!");
	fflush(stdout);

	tcoords[0]=coords[0]+1;
	MPI_Cart_rank(torus,tcoords,&trank);
	MPI_Isend(lay[bS+1],bS+2,MPI_CHAR,trank,8,torus,&request);
	MPI_Recv(lay[0],bS+2,MPI_CHAR,MPI_ANY_SOURCE,8,torus,&status);
	MPI_Wait(&request,&status);
	printf("top gostsent!");
	fflush(stdout);
}
void makeStep(char **lay1,char **lay2,char *lc,char *rc,int bS){
	char sum;
	for(int i=1;i<bS+1;++i){
		for(int j=1;j<bS+1;++j){
			sum=lay1[i-1][j-1]+lay1[i-1][j]+lay1[i-1][j+1]+lay1[i][j-1]+lay1[i][j+1]+lay1[i+1][j-1]+lay1[i+1][j]+lay1[i+1][j+1];
			if(sum==3)lay2[i][j]=1;
			else if(sum>3 || sum<2)lay2[i][j]=0;
		}
		lc[i-1]=lay1[i][1];
		rc[i-1]=lay1[i][bS];
	}
}
int logChanges(char **lay1,char **lay2,int bS,int coords[2],int iter){
	int res=0;
	for(int i=1;i<bS+1;++i){
		for(int j=1;j<bS+1;++j){
			if(lay1[i][j]!=lay2[i][j]){
				res=1;
				printf("i=%d p=%d;%d c=%d;%d\n",iter,coords[0],coords[1],i-1,j-1);
			}
		}
	}
	return res;
}
int main(int argc, char *argv[]){
	int size, rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm torus;
	int matrixSize=atoi(argv[1]);
	int blockSize,cartSize;//blockSize - размер блока матрицы cartSize - размер декартовой топологии
	//Проверяем размерность топологии и блоков матрицы (обе должны быть квадратными)
	if(isInteger(sqrt((double)size),cartSize) && isInteger(sqrt(((double)matrixSize*matrixSize)/size),blockSize)){
		if(rank==0)printf("Block size=%d\nCart size=%d\n",blockSize,cartSize);
		//Init
		MPI_Request request;
		MPI_Status status;
		int coords[2];
		int tcoords[]={0,0};
		int trank;
		{
			int dims[]={cartSize,cartSize};
			int periods[]={true,true};
			MPI_Cart_create(MPI_COMM_WORLD,2,dims,periods,true,&torus);
			int nrank;
			MPI_Comm_rank(torus,&nrank);
			MPI_Cart_coords(torus,nrank,2,coords);
		}
		char **lay1,**lay2,*leftC,*rightC;
		leftC=new char [blockSize];
		rightC=new char [blockSize];
		memset(leftC,0,blockSize);
		memset(rightC,0,blockSize);
		lay1=new char* [blockSize+2];
		lay2=new char* [blockSize+2];
		for(int i=0;i<blockSize+2;++i){
			lay1[i]=new char[blockSize+2];
			memset(lay1,0,blockSize+2);
			printf("%d",lay1[0][0]);
			lay2[i]=new char[blockSize+2];
			memset(lay2,0,blockSize+2);
		}
		
		char *block=new char[blockSize*blockSize];
		if(coords[0]==0 && coords[1]==0){
			generateBlock(block,blockSize);
			printf("\np=%d;%d\n",coords[0],coords[1]);
			printBlock(block,blockSize);
			blockToLay(block,lay1,leftC,rightC,blockSize);
			tcoords[1]=1;
			for(tcoords[0]=0;tcoords[0]<cartSize;++tcoords[0]){
				for(tcoords[1];tcoords[1]<cartSize;++tcoords[1]){
					generateBlock(block,blockSize);
					printf("\np=%d;%d\n",tcoords[0],tcoords[1]);
					printBlock(block,blockSize);
					MPI_Cart_rank(torus,tcoords,&trank);
					MPI_Send(block,blockSize*blockSize,MPI_CHAR,trank,1,torus);
				}
				tcoords[1]=0;
			}
		}else{
			MPI_Recv(block,blockSize*blockSize,MPI_CHAR,MPI_ANY_SOURCE,1,torus,&status);
			blockToLay(block,lay1,leftC,rightC,blockSize);
		}
		delete []block;
		ghostShare(torus,lay1,leftC,rightC,coords,blockSize,cartSize);
		printf("\nI=%d;%d\n",coords[0],coords[1]);
		fflush(stdout);
		MPI_Barrier(torus);
		//Start
		long iter;
		int changes;
		for(iter=0;iter<10000000;++iter){
			if(iter%2==0){
				makeStep(lay1,lay2,leftC,rightC,blockSize);
			}else{
				makeStep(lay2,lay1,leftC,rightC,blockSize);
			}
			changes=logChanges(lay1,lay2,blockSize,coords,iter);
			MPI_Allreduce(&changes,&changes,1,MPI_INTEGER,MPI_MAX,torus);
			if(changes==0)break;
			else{
				if(iter%2==0)ghostShare(torus,lay2,leftC,rightC,coords,blockSize,cartSize);
				else ghostShare(torus,lay1,leftC,rightC,coords,blockSize,cartSize);
			}
		}
		delete[] leftC;
		delete[] rightC;
		for(int i=0;i<blockSize+2;++i){
			delete[] lay1[i];
			delete[] lay2[i];
		}
		delete[] lay1;
		delete[] lay2;
	}else if(rank==0){
		printf("Size is incorrect!");
	}
	MPI_Finalize();
	return 0;
}