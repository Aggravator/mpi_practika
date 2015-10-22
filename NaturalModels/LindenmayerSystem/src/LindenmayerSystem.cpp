/*
 ============================================================================
 Name        : LindenmayerSystem.c
 Author      : 
 Version     :
 Copyright   : Your copyright notice
 Description : Compute Pi in MPI C++
 ============================================================================
 */
#include <mpi.h>
#include <cstdlib>
#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>

using namespace std;
typedef multimap<char, pair<double,string> > RuleMap;
string getValue(char key,RuleMap &R,double p){
	pair<RuleMap::iterator,RuleMap::iterator> range=R.equal_range(key);
	if(range.first==range.second)return string()+key;
	RuleMap::iterator i=range.first;
	double pSum=0;

	do{
		pSum+=i->second.first;
	}while(pSum<p && ++i!=range.second);
	if(pSum>p) return i->second.second;
	else return string()+key;

}
char* exchange(string& buffer,int from,int to){
	MPI_Status status;
	MPI_Request request;
	unsigned length=buffer.length(),neighbourLength;
	MPI_Sendrecv(&length,1,MPI_INT,from,1,&neighbourLength,1,MPI_INT,to,1,MPI_COMM_WORLD,&status);
	int c=length>neighbourLength ? (length-neighbourLength)/2 : 0;
	MPI_Isend(buffer.data()+length-c,c,MPI_CHAR,to,2,MPI_COMM_WORLD,&request);

	char *temp;
	if(from!=MPI_PROC_NULL){
		MPI_Probe(to,2,MPI_COMM_WORLD,&status);
		int count=0;
		MPI_Get_count(&status,MPI_CHAR,&count);
		char *temp=new char[count +1];
		temp[count]='\0';
		MPI_Recv(temp,count,MPI_CHAR,to,2,MPI_COMM_WORLD,&status);
	}else{
		temp=new char[1];
		temp[0]='\0';
	}

	MPI_Wait(&request,&status);
	return temp;
}
void align(string& buffer,int rank,int size){
	char *temp;
	temp=exchange(buffer,rank!=0 ?(rank-1)%size:MPI_PROC_NULL,rank!=size-1 ?(rank+1)%size:MPI_PROC_NULL);
	buffer.insert(0,temp);
	delete[] temp;

	temp=exchange(buffer,rank!=size-1 ?(rank+1)%size:MPI_PROC_NULL,rank!=0 ?(rank-1)%size:MPI_PROC_NULL);
	buffer+=temp;
	delete[] temp;
}

int main(int argc, char *argv[]) {
	int rank, size;
	MPI_Status status;
	double time;
	vector<unsigned> statistic;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);

	int T=4,k=20;
	if(argc>1){
		T = atoi(argv[1]);
		k=atoi(argv[2]);
	}

	RuleMap R;
	R.insert(make_pair('a',make_pair(1.0,"b")));
	R.insert(make_pair('b',make_pair(1.0,"ab")));

	int bufferCount=2;
	string *buffers=new string[bufferCount];
	buffers[0]="a";
	int curBuffer=0,nextBuffer;

	if(rank==0) time=MPI_Wtime();

	for(int t=0; t<T; t++ ){
		nextBuffer=(curBuffer+1)%bufferCount;
		buffers[nextBuffer].clear();
		for(unsigned i=0;i<buffers[curBuffer].length();++i){
			buffers[nextBuffer]+=getValue(buffers[curBuffer][i],R,double(rand())/RAND_MAX);
		}
		curBuffer=nextBuffer;
		if(t%k==0 && t!=0){
			statistic.push_back(buffers[curBuffer].length());
			align(buffers[curBuffer],rank,size);
		}
	}

	unsigned *lengths=new unsigned[size];
	unsigned length=buffers[curBuffer].length();
	MPI_Allgather(&length,1,MPI_INT,lengths,1,MPI_INT,MPI_COMM_WORLD);
	int pos=0;
	for(int i=0;i<rank;++i)pos+=lengths[i];
	MPI_File f;
	MPI_File_open(MPI_COMM_WORLD, "output.txt", MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &f);
	MPI_File_set_view(f, pos, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);
	MPI_File_write(f, buffers[curBuffer].data(), buffers[curBuffer].length(), MPI_CHAR, &status);
	MPI_File_close(&f);
	delete[] lengths;

	delete[] buffers;

	if(rank==0){
		unsigned* statisticInf=new unsigned[statistic.size()*size];
		MPI_Gather(statistic.data(),statistic.size(),MPI_INT,statisticInf,statistic.size(),MPI_INT,rank,MPI_COMM_WORLD);

		unsigned sum;
		ofstream f("stat.txt");
		f<<"      ";
		for(int i=0;i<size;++i)f<<setw(10)<<i;
		f<<endl;
		for(int i=0;i<statistic.size();++i){
			f<<setw(6)<<(i+1)*k;
			sum=0;
			for(int j=0;j<size;++j)sum+=statisticInf[j*statistic.size()+i];
			for(int j=0;j<size;++j)
				f<<setw(10)<<statisticInf[j*statistic.size()+i]/double(sum);
			f<<endl;
		}
		f.close();
		delete[] statisticInf;
	}else{
		MPI_Gather(statistic.data(),statistic.size(),MPI_INT,void,statistic.size(),MPI_INT,0,MPI_COMM_WORLD);
	}

	MPI::Finalize();
	return 0;
}

