#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <vector>
#include <omp.h>

#define RAND_MAX 6074
static int next = 1;

//double abs(double value){return value>0 ? value : -value;}

unsigned randm(void)
{
	const static int a = 106, c = 1283, m = 6075;
	next = (a*next + c) % m;
	return next;
}
void srandm(unsigned v){
	next = v%RAND_MAX;
}
int getPeriod(unsigned(*s) (void)){
	std::vector<int>elements;
	int i, j;
	for (i = 0; i<6077; ++i){
		elements.push_back(s());
		if (elements[i] == elements[0] && i != 0)return i;
	}
	return i;
}
double parabolaFunc(double x){
	return -x*x + 5 * x + 2;
}
double integralPF(double x){
	return (x*(x*(-2 *x + 15) + 12)) / 6.0;

}
double sequential(double x1, double x2, double h, double(*func)(double x), double eps, int &iter){
	int areas = 0, iteration = 0;
	double computedArea = 0,oldComputedArea=0;
	do{
		++iteration;
		if (func(x1 + randm()*(x2-x1) / RAND_MAX)>(randm()*h / RAND_MAX))++areas;
		oldComputedArea = computedArea;
		computedArea = (x2 - x1)*h*double(areas) / double(iteration);
	} while (abs(oldComputedArea-computedArea)>eps || iteration<20);
 	iter = iteration;
	return computedArea;
}
double parallel(double x1, double x2, double h, double(*func)(double x), double eps, int &iter, unsigned short thCount){
	double computedArea = 0;
	int iteration=0;
#pragma omp parallel num_threads(thCount) reduction(+:computedArea,iteration)
	{
		int areas = 0;
		iteration = 0;
		double oldComputedArea = 0;
		double length = (x2 - x1) / omp_get_num_threads();
		double startx = length*omp_get_thread_num() + x1;
		double realArea = integralPF(startx + length) - integralPF(startx);
		computedArea = 0;
		do{
			++iteration;
			if (func(startx + randm()*(length) / RAND_MAX)>(randm()*h / RAND_MAX))++areas;
			oldComputedArea = computedArea;
			computedArea = (length)*h*areas / double(iteration);
		} while (abs(oldComputedArea - computedArea)>eps || iteration<20);
	}
	iter = iteration / thCount;
	return computedArea;
}
int main(){
	srandm(time(0));
	int period=getPeriod(randm);
	printf("Period: %d\n",period);
	int k=period;
	unsigned expValue=period/k;
	std::vector<int> hiarr(k);
	std::fill(hiarr.begin(),hiarr.end(),0);
	double elementink=static_cast<double>(period)/k;
	for(int i=0;i<period;++i) hiarr[int(randm()/elementink)]+=1;
	double hi2=0;
	for(int i=0;i<k;++i)hi2+=(hiarr[i]-expValue)*(hiarr[i]-expValue)/expValue;
	printf("hi2: %lf\n",hi2);

	double x1, x2, h;
	//-x^2+5x+2
	x1 = (-5 + sqrt(33.0)) / -2.0;
	x2 = (-5 - sqrt(33.0)) / -2.0;
	h = 8.5;
	//h = 8.25;

	double area = (x2 - x1)*h;
	double parabolaArea = integralPF(x2) - integralPF(x1);
	printf("\nreal area: %lf\n", parabolaArea);
	double computeS, time,sumTime=0,minTime=10000,maxTime=-10;
	int iter,i;
	unsigned long sumIter = 0;
	const int repCount = 15;
	for (i = 0; i < repCount; ++i){
		time = omp_get_wtime();
		computeS = sequential(x1, x2, h, parabolaFunc, 0.000001, iter);
		time = omp_get_wtime() - time;
		if (time < minTime)minTime = time;
		if (time > maxTime)maxTime = time;
		sumTime += time;
		sumIter += iter;
	}
	printf("\nthread count: 1\n");
	printf("time: %lf\nminTime:%lf\tmaxTime:%lf\n", sumTime/i,minTime,maxTime);
	printf("computed area: %lf\n", computeS);
	printf("iteration count: %d\n", sumIter/i);
	fflush(stdout);
	for (int j = 2; j <= 16; ++j){
		sumIter = 0;
		sumTime = 0;
		minTime = 10000; maxTime = -10;
		for (i = 0; i < repCount; ++i){
			time = omp_get_wtime();
			computeS = parallel(x1, x2, h, parabolaFunc, 0.000001, iter, j);
			time = omp_get_wtime() - time;
			if (time < minTime)minTime = time;
			if (time > maxTime)maxTime = time;
			sumTime += time;
			sumIter += iter;
		}
		printf("\nthread count: %d\n", j);
		printf("time: %lf\nminTime:%lf\tmaxTime:%lf\n", sumTime / i, minTime, maxTime);
		printf("computed area: %lf\n", computeS);
		printf("iteration count: %d\n", sumIter / i);
		fflush(stdout);
	}
	return 0;
}