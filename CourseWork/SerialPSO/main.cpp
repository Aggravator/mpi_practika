#include <iostream>
#include <cstdlib>
#include <vector>
#include <limits>
#include <conio.h>
#include <algorithm>
#include <iomanip>
#include <chrono>
using namespace std;
const double PI = 3.1415926535897932384626433832795;
inline double rndDouble(double min, double max){
	return (rand() / double(RAND_MAX))*(max - min) + min;
}
const int dimensionCount = 10;
struct Particle{
	double fitness;
	double bestFitness;
	double pos[dimensionCount];
	double velocity[dimensionCount];
	double bestPos[dimensionCount];
};
double fitnessFunc(const double *pos){
	//Rastrigin's function
	/*double temp = 0;
	for (int i = 0; i < dimensionCount; ++i)
		temp += (pos[i] * pos[i] - 10 * cos(2 * PI*pos[i]));
	return 10 * dimensionCount + temp;*/
	//Schwefel's function
	double temp = 0;
	for (int i = 0; i < dimensionCount; ++i)temp += (-pos[i] * sin(sqrt(abs(pos[i]))));
	return temp;
}
double fitnessFunc(const Particle &particle){
	return fitnessFunc(particle.pos);
}
double posCmp(const double *pos1, const double *pos2){
	double temp = 0;
	for (int i = 0; i < dimensionCount; ++i)temp += abs(pos2[i]-pos1[i]);
	return temp;
}
bool isEqual(double a, double b, double delta){
	if (abs(a - b) < delta)return true;
	else return false;
}

int main(){
	srand(time(0));
	int particleCount = 500;
	vector<Particle> particles(particleCount);
	particles.reserve(particleCount);
	double minDimValue=-500, maxDimValue=500;
	double c1=2, c2=2;//velocity coefficients personal and global
	double w =0.8,wMax=1,wMin=0.7;
	double bestFitness,oldBestFitness;
	double bestPos[dimensionCount], oldBestPos[dimensionCount];
	double *localBestPos;
	double delta = 0.001;
	double deltaForCdiw=0.0001;
	double z;
	int repeatCounter;
	int bestPosIndex;
	int testCount = 0,rightCount=0;
	int iteration;
	int maxIterations = 1200, maxRepeatedIter = 500;
	chrono::high_resolution_clock::time_point t1, t2;
	unsigned int sumTime = 0;
again:	
	t1 = chrono::high_resolution_clock::now();
	//Initialize 
	repeatCounter = 0;
	/*do{
		z = rndDouble(0.001, 0.999);
	} while (isEqual(z, 0.25, deltaForCdiw) || isEqual(z, 0.5, deltaForCdiw) || isEqual(z, 0.75, deltaForCdiw));*/

	oldBestFitness=bestFitness = numeric_limits<double>::max();
	for (int i = 0; i < particleCount; ++i){
		for (int j = 0; j < dimensionCount; ++j){
			particles[i].pos[j] = rndDouble(minDimValue, maxDimValue);
			//particles[i].velocity[j] = rndDouble(0, (maxDimValue - minDimValue) / 6.0);
		}
		copy(particles[i].pos, particles[i].pos + dimensionCount, particles[i].bestPos);
		fill(particles[i].velocity, particles[i].velocity + dimensionCount, 0);
		particles[i].fitness = fitnessFunc(particles[i]);
		particles[i].bestFitness = particles[i].fitness;
		if (bestFitness>particles[i].fitness){
			bestFitness = particles[i].fitness;
			copy(particles[i].pos, particles[i].pos + dimensionCount,bestPos);
		}
	}
	copy(bestPos, bestPos + dimensionCount, oldBestPos);
	//main cycle
	iteration = 0;
	bool ju = true;
	while (/*iteration < maxIterations &&*/ repeatCounter<maxRepeatedIter || iteration<100*dimensionCount){
		++iteration;
		bestPosIndex = -1;
		//w = wMax - iteration*(wMax - wMin) / double(maxIterations);
		w = (wMax - wMin)*(maxIterations - iteration) / maxIterations + wMin;
		//z = 4 * z*(1 - z);w = (wMax - wMin)*(maxIterations - iteration) / maxIterations + wMin*z;
		for (int i = 0; i < particleCount; ++i){
			localBestPos=min({ particles[i], particles[(i + particleCount - 1) % particleCount], particles[(i + particleCount + 1) % particleCount] },
				[](const Particle& a,const Particle& b){return a.bestFitness < b.bestFitness; }).bestPos;
			//ñhange of velocity and position of particle
			for (int j = 0; j < dimensionCount; ++j){
				//particles[i].velocity[j] = w*particles[i].velocity[j] + c1*rndDouble(0, 1)*(particles[i].bestPos[j] - particles[i].pos[j]) + c2*rndDouble(0, 1)*(bestPos[j] - particles[i].pos[j]);
				particles[i].velocity[j] = w*particles[i].velocity[j] + c1*rndDouble(0, 1)*(particles[i].bestPos[j] - particles[i].pos[j]) + c2*rndDouble(0, 1)*(localBestPos[j] - particles[i].pos[j]);
				//particles[i].velocity[j] = 0.729*(particles[i].velocity[j] + 2.05*rndDouble(0, 1)*(particles[i].bestPos[j] - particles[i].pos[j]) + 2.05*rndDouble(0, 1)*(localBestPos[j] - particles[i].pos[j]));
				//particles[i].velocity[j] = 0.729*(particles[i].velocity[j] + 2.05*rndDouble(0, 1)*(particles[i].bestPos[j] - particles[i].pos[j]) + 2.05*rndDouble(0, 1)*(bestPos[j] - particles[i].pos[j]));
				particles[i].pos[j] += particles[i].velocity[j];
				//particle position correction
				if (particles[i].pos[j] < minDimValue){
					particles[i].pos[j] = minDimValue;
					particles[i].velocity[j] = 0;
				}
				if (particles[i].pos[j] > maxDimValue){
					particles[i].pos[j] = maxDimValue;
					particles[i].velocity[j] = 0;
				}
			}
			double tfitness = fitnessFunc(particles[i]);
			//correction the best position of particle
			if (particles[i].bestFitness > tfitness){
				copy(particles[i].pos, particles[i].pos + dimensionCount, particles[i].bestPos);
				particles[i].bestFitness = tfitness;
			}
			particles[i].fitness = tfitness;
			//search new the best search position
			if (tfitness < bestFitness){
				bestPosIndex = i;
				bestFitness = tfitness;
			}
		}
		//correction the best global position
		//if (bestPosIndex != -1) copy(particles[bestPosIndex].pos, particles[bestPosIndex].pos + dimensionCount, bestPos);
		//auto worseE = min_element(particles.begin(), particles.end(), [](const Particle& p1, const Particle& p2){return p1.fitness>p2.fitness; });
		//cout <<iteration<< ": Best|Worst fitness: " <<setw(13)<< fitnessFunc(bestPos) <<" | "<<setw(13)<<worseE->fitness<< endl;
		if (/*posCmp(bestPos,oldBestPos)>delta*/abs(bestFitness-oldBestFitness)>delta){
			repeatCounter = 0;
			//copy(bestPos, bestPos + dimensionCount, oldBestPos);
			oldBestFitness = bestFitness;
		}
		else ++repeatCounter;
		if (maxIterations < iteration + 2)maxIterations += 40;
	}
	t2 = chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count()/1000000.0;
	sumTime += std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
	cout << "Launch #: " << testCount << endl;
	cout << "Best position: [ "<<setprecision(5)<<bestPos[0];
	for (int i = 1; i < dimensionCount; ++i)cout << ", " << bestPos[i];
	cout << "]" << endl;
	cout << "Best fitness: " << bestFitness << endl;
	cout << "Time: " << duration<<endl;
	cout << "Iterations: " << iteration << endl;
	cout << "Avg answers: " << double(rightCount) / testCount << endl;
	cout << "---------------------------------------" << endl;
	if (abs(bestFitness + dimensionCount*418.9829) < 0.1)rightCount += 1;
	++testCount;
	if (testCount<500) goto again;
	cout << "Launch count: " << testCount<<endl;
	cout << "Right answers: " << rightCount <<" | "<<double(rightCount)/testCount<< endl;
	cout << "Average time: " << (sumTime / double(testCount)) / 1000000.0 << endl;
	getch();
}