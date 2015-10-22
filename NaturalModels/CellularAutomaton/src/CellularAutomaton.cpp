#include <mpi.h>
#include <iostream>
#include <cstdlib>
#include <bitset>
#include <fstream>

using namespace std;

int main(int argc, char *argv[]) {
        int rank, size;
        MPI_Status status;
        double time;

        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);
        MPI_Comm_size(MPI_COMM_WORLD,&size);

        int alpha=90;
        int n=19;
        int T=3;

        if(argc>1){
                alpha = atoi(argv[1]);
                n = atoi(argv[2]);
                T = atoi(argv[3]);

        }

        bitset<8>  R;
        for( int i=0; i<8; i++ )
        	R[i] = alpha & (1<<i);

        char *data=new char[n+2];
        fill(data,data+n+2,false);

        if(rank==0) time=MPI_Wtime();

        int mid=(n*size+2)/2;
        if(mid>n*rank && mid<=n*(rank+1)+1) data[mid-n*rank]=1;

        char a=data[0],b=data[1];
        for(int t=0;t<T;++t){
                for(int j=1;j<n+1;++j){
                        data[j]=R[(a<<2)+(b<<1)+data[j+1]];
                        a=b;
                        b=data[j+1];
                }

                MPI_Sendrecv(data+1,1,MPI_CHAR,rank-1>=0 ? rank-1 : MPI_PROC_NULL,1,data+n+1,1,MPI_CHAR,rank+1<size ? rank+1 : MPI_PROC_NULL,1,MPI_COMM_WORLD,&status);
                MPI_Sendrecv(data+n,1,MPI_CHAR,rank+1<size ? rank+1 : MPI_PROC_NULL,2,data,1,MPI_CHAR,rank-1>=0 ? rank-1 : MPI_PROC_NULL,2,MPI_COMM_WORLD,&status);
        }

        for(int i=0;i<n+2;++i)data[i]+='1'-1;

        MPI_File f;
        MPI_File_open(MPI_COMM_WORLD, "output.txt", MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &f);
        MPI_File_set_view(f, rank*n, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);
        MPI_File_write(f, data+1, n, MPI_CHAR, &status);
        MPI_File_close(&f);

        if(rank==0){
                ofstream f("stat.txt");
                f<<"Time: "<<MPI_Wtime()-time<<endl;
                f<<"Proc count: "<<size<<endl;
                f<<"Alpha: "<<alpha<<endl;
                f<<"Iterations: "<<T<<endl;
                f<<"n(N): "<<n<<"("<<n*size+2<<")"<<endl;
                f.close();
        }

        delete[] data;
        MPI::Finalize();
        return 0;
}


