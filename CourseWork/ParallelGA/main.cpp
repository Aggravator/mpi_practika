
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <omp.h>


using namespace std;

inline double frand(double a, double b)
{
	return a + (b - a)*(rand() / double(RAND_MAX));
}
inline bool isFileExist(const char* name) {
	ifstream f(name);
	if (f.good()) {
		f.close();
		return true;
	}
	else {
		f.close();
		return false;
	}
}

inline int do_walk(int a, int b, int x, double p, int& t)
{
	int step = 0;
	while (x>a && x<b)
	{
		if (frand(0, 1)<p)
			x += 1;
		else
			x -= 1;
		t += 1.0;
		step += 1;
	}
	return x;
}


int main(int argc, char** argv)
{
	int a = -100, b = 100, x = 0, N = 20000;
	int threadCount = 1;
	double p = 0.5;

	if (argc >= 6){
		a = atoi(argv[1]);
		b = atoi(argv[2]);
		x = atoi(argv[3]);
		p = atof(argv[4]);
		N = atoi(argv[5]);
		threadCount = atoi(argv[6]);
	}

	int t = 0;
	int w = 0;

	//omp_set_num_threads(threadCount);
	double time = omp_get_wtime();
#pragma omp parallel for schedule(dynamic,100) reduction(+:w,t) firstprivate(a,b,x,p,N)  if(threadCount>1) num_threads(threadCount)
	for (int i = 0; i<N; i++)
	{
		int cx = x;
		while (cx>a && cx<b)
		{
			if (frand(0, 1)<p)
				cx += 1;
			else
				cx -= 1;
			t += 1;
		}
		if (cx == b)
			w += 1;
	}
	time = omp_get_wtime() - time;

	string filename("output");
	int fileN = 1;
	while (isFileExist((filename + to_string(fileN) + ".txt").c_str()))++fileN;

	ofstream f(filename + to_string(fileN) + ".txt");
	f << double(w) / N << " " << double(t) / N << endl;
	f.close();

	f.open(string("stat") + to_string(fileN) + ".txt");
	f << "Time: " << time << endl;
	f << "Thread count: " << threadCount << endl;
	f << "a: " << a << endl;
	f << "b: " << b << endl;
	f << "x: " << x << endl;
	f << "p: " << p << endl;
	f << "N: " << N << endl;
	f << "Time: " << time << endl;
	cout << "Time: " << time << endl;
	cout << "Thread count: " << threadCount << endl;
	cout << double(w) / N << " " << double(t) / N << endl;
	f.close();
	cin >> filename;
	return 0;
}

