#include <iostream>
#include <mpi.h>
#include <cstdlib>
#include <vector>
#include <limits>
#include <algorithm>
#include <iomanip>
#include <cmath>
using std::vector;
using std::cout;
const double PI = 3.1415926535897932384626433832795;
double dRand(double fMin, double fMax){
	return fMin + (double(rand()) / RAND_MAX) * (fMax - fMin);
}
/*bool isEqual(double a, double b, double delta){
	if (std::abs(a - b) < delta)return true;
	else return false;
}*/
int rank, size;
struct {
	int particleCount;
	int dimCount;
}params;
struct Particle{
	double fitness;
	double bestFitness;
	double *pos;
	double *velocity;
	double *bestPos;
	Particle(){
		pos = new double[params.dimCount];
		bestPos = new double[params.dimCount];
		velocity = new double[params.dimCount];
	}
	Particle(const Particle &p){
		pos = new double[params.dimCount];
		bestPos = new double[params.dimCount];
		velocity = new double[params.dimCount];
		std::copy(p.pos, p.pos + params.dimCount, pos);
		std::copy(p.bestPos, p.bestPos + params.dimCount, bestPos);
		std::copy(p.velocity, p.velocity + params.dimCount, velocity);
		fitness = p.fitness;
		bestFitness = p.bestFitness;
	}
	Particle& operator=(const Particle& p){
		std::copy(p.pos, p.pos + params.dimCount, pos);
		std::copy(p.bestPos, p.bestPos + params.dimCount, bestPos);
		std::copy(p.velocity, p.velocity + params.dimCount, velocity);
		fitness = p.fitness;
		bestFitness = p.bestFitness;
		return *this;
	}
	void indToCh(char* dest){
		memcpy(dest, &fitness, sizeof(fitness));
		dest += sizeof(fitness);
		memcpy(dest, pos, sizeof(double)*params.dimCount);
	}
	void chToInd(const char* dest){
		memcpy(&fitness, dest, sizeof(fitness));
		dest += sizeof(fitness);
		memcpy(pos, dest, sizeof(pos[0])*params.dimCount);
	}
};
bool bestFitParticle(const Particle& p1, const Particle& p2){
	return p1.bestFitness < p2.bestFitness;
}
double rastriginF(const double *pos){//Rastrigin's function
	double temp = 0;
	for (int i = 0; i < params.dimCount; ++i)
		temp += (pos[i] * pos[i] - 10 * cos(2 * PI*pos[i]));
	return 10 * params.dimCount + temp;
}
bool rastriginFCheck(const double answer, const double delta = 0.001){
	return std::abs(answer - 0) < delta;
}
double schwefelF(const double *pos){//Schwefel's function
	double temp = 0;
	for (int i = 0; i < params.dimCount; ++i)temp += (-pos[i] * sin(sqrt(std::abs(pos[i]))));
	return temp;
}
bool scwefelFCheck(const double answer, const double delta = 0.001){
	return std::abs(answer - params.dimCount*418.9829) < delta;
}
double rosenbrockF(const double *pos){//Rosenbrock's function
	double temp = 0;
	for (int i = 0; i < params.dimCount - 1; ++i){
		temp += 100 * (pos[i + 1] - pos[i] * pos[i])*(pos[i + 1] - pos[i] * pos[i]) + (pos[i] - 1)*(pos[i] - 1);
	}
	return temp;
}
bool rosenbrockFCheck(const double answer, const double delta = 0.001){
	return std::abs(answer - 0) < delta;
}
double griewangkF(const double *pos){
	double temp1 = 0, temp2 = 1;
	for (int i = 0; i < params.dimCount; ++i){
		temp1 += pos[i] * pos[i] / 4000.0;
		temp2 *= cos(pos[i] / sqrt(double(i+1)));
	}
	return 1 + temp1 - temp2;
}
bool griewangkFCheck(const double answer, const double delta = 0.001){
	return std::abs(answer - 0) < delta;
}
class Task{
public:
	double minValue, maxValue;
	double(*fitness)(const double *genotype);
	bool(*check)(const double answer, const double delta);
};

void taskInit(Task& task, int taskId){
	if (taskId == 0){//rosenbrock
		task.fitness = rosenbrockF;
		task.check = rosenbrockFCheck;
		task.minValue = -2.048;
		task.maxValue = 2.048;
	}
	else if (taskId == 1){//rastrigin
		task.fitness = rastriginF;
		task.check = rastriginFCheck;
		task.minValue = -5.12;
		task.maxValue = 5.12;
	}
	else if (taskId == 2){//schwefel
		task.fitness = schwefelF;
		task.check = scwefelFCheck;
		task.minValue = -500;
		task.maxValue = 500;
	}
	else if (taskId == 3){//griewangk
		task.fitness = griewangkF;
		task.check = griewangkFCheck;
		task.minValue = -600;
		task.maxValue = 600;
	}
}
char tasksName[4][15] = { "Rosenbrock", "Rastrigin", "Schwefel", "Griewangk" };

int main(int argc, char **argv){
	int testCount = 100, successTests = 0, testIter, maxIter = 2000, iteration, taskId = 0, maxRepeatIter = 200, repeatIter;
	double delta = 0.001;
	double t1, t2;
	//params initializing
	params.particleCount = 100;
	params.dimCount = 10;
	{
		int tempi;
		double tempd;
		bool isElitCountCh = false;
		for (int i = 0; i < argc; ++i){
			if (strncmp(argv[i], "pc=", strlen("pc=")) == 0){
				tempi = atoi(argv[i] + strlen("pc="));
				if (tempi != 0)	params.particleCount = tempi;
			}
			else if (strncmp(argv[i], "dc=", strlen("dc=")) == 0){
				tempi = atoi(argv[i] + strlen("dc="));
				if (tempi != 0) params.dimCount = tempi;
			}
			else if (strncmp(argv[i], "tc=", strlen("tc=")) == 0){
				tempi = atoi(argv[i] + strlen("tc="));
				if (tempi != 0)	testCount = tempi;
			}		
			else if (strncmp(argv[i], "mi=", strlen("mi=")) == 0){
				tempi = atof(argv[i] + strlen("mi="));
				if (tempi != 0) maxIter = tempi;
			}
			else if (strncmp(argv[i], "ta=", strlen("ta=")) == 0){
				tempi = atoi(argv[i] + strlen("ta="));
				if (tempi != 0)	taskId = tempi;
			}
		}
	}
	srand(time(0));
	vector<Particle> particles(params.particleCount);
	bool isTerminate;
	double oldFitness, bestFitness;
	Task task;
	taskInit(task, taskId);

	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0){
		cout << "Parameters:" << std::endl;
		cout << "\tParticle count: " << params.particleCount << std::endl;
		cout << "\tDimension count: " << params.dimCount << std::endl;
		cout << "\tTask: " << tasksName[taskId] << std::endl;
		cout << "\tNode count:" << size << std::endl;
	}
	for (testIter = 0; testIter < testCount; ++testIter){
		if (rank == 0)t1 = MPI_Wtime();
		bestFitness = std::numeric_limits<double>::max();
		iteration = repeatIter=0;
		isTerminate = false;
		//Initialize 
		for (int i = 0; i < params.particleCount; ++i){
			for (int j = 0; j < params.dimCount; ++j)
				particles[i].pos[j] = dRand(task.minValue, task.maxValue);
			std::copy(particles[i].pos, particles[i].pos + params.dimCount, particles[i].bestPos);
			std::fill(particles[i].velocity, particles[i].velocity + params.dimCount, 0);
			particles[i].fitness = task.fitness(particles[i].pos);
			particles[i].bestFitness = particles[i].fitness;
			if (bestFitness>particles[i].fitness)bestFitness = particles[i].fitness;
		}
		oldFitness = bestFitness;
		//main cycle
		while (true){
			++iteration;
			for (int i = 0; i < params.particleCount; ++i){
				const double *localBestPos = std::min(std::min(particles[i], particles[(i + params.particleCount - 1) % params.particleCount], bestFitParticle), particles[(i + params.particleCount + 1) % params.particleCount], bestFitParticle).bestPos;
				//ñhange of velocity and position of particle
				for (int j = 0; j < params.dimCount; ++j){
					particles[i].velocity[j] = 0.729*(particles[i].velocity[j] + 2.05*dRand(0, 1)*(particles[i].bestPos[j] - particles[i].pos[j]) + 2.05*dRand(0, 1)*(localBestPos[j] - particles[i].pos[j]));
					particles[i].pos[j] += particles[i].velocity[j];
					//particle position correction
					if (particles[i].pos[j] < task.minValue){
						particles[i].pos[j] = task.minValue;
						particles[i].velocity[j] = 0;
					}
					else if (particles[i].pos[j] > task.maxValue){
						particles[i].pos[j] = task.maxValue;
						particles[i].velocity[j] = 0;
					}
				}
				double tfitness = task.fitness(particles[i].pos);
				//correction the best position of particle
				if (particles[i].bestFitness > tfitness){
					if (bestFitness > tfitness)bestFitness = tfitness;
					std::copy(particles[i].pos, particles[i].pos + params.dimCount, particles[i].bestPos);
					particles[i].bestFitness = tfitness;
				}
				particles[i].fitness = tfitness;
			}
			if (rank == 0 && std::abs(oldFitness - bestFitness)<delta)++repeatIter;
			else{
				repeatIter = 0;
				oldFitness = bestFitness;
			}
			if (iteration % 10 == 0){
				if (rank == 0){
					isTerminate = repeatIter>maxRepeatIter-1 || iteration > maxIter-1;
					MPI_Bcast(&isTerminate, sizeof(bool), MPI_CHAR, 0, MPI_COMM_WORLD);
				}
				else{
					MPI_Bcast(&isTerminate, sizeof(bool), MPI_CHAR, 0, MPI_COMM_WORLD);
				}
				if (isTerminate)break;
				Particle *best = &particles[0], *worst1 = &particles[0], *worst2 = &particles[0];
				for (int i = 1; i < params.particleCount; ++i){
					if (particles[i].bestFitness>worst1->bestFitness){
						worst2 = worst1;
						worst1 = &particles[i];
					}
					if (particles[i].bestFitness < best->bestFitness)best = &particles[i];
				}
				if (size > 1){
					if (rank % 2 == 0){
						MPI_Sendrecv(best->pos, params.dimCount, MPI_DOUBLE, (rank + 1) % size, 0, worst1->pos, params.dimCount, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
						MPI_Sendrecv(best->pos, params.dimCount, MPI_DOUBLE, (rank + size - 1) % size, 0, worst2->pos, params.dimCount, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
					}
					else{
						MPI_Sendrecv(best->pos, params.dimCount, MPI_DOUBLE, (rank + size-1) % size, 0, worst1->pos, params.dimCount, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
						MPI_Sendrecv(best->pos, params.dimCount, MPI_DOUBLE, (rank + 1) % size, 0, worst2->pos, params.dimCount, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
					}
					worst1->fitness = task.fitness(worst1->pos);
					worst2->fitness = task.fitness(worst2->pos);
					if (worst1->fitness < worst1->bestFitness){
						std::copy(worst1->pos, worst1->pos + params.dimCount, worst1->bestPos);
						worst1->bestFitness = worst1->fitness;
					}
					if (worst2->fitness < worst2->bestFitness){
						std::copy(worst2->pos, worst2->pos + params.dimCount, worst2->bestPos);
						worst2->bestFitness = worst2->fitness;
					}
				}
			}
		}
		{
			Particle *best = &particles[0];
			for (int i = 1; i < params.particleCount; ++i)
				if (particles[i].bestFitness < best->bestFitness)best = &particles[i];
			int sizeI = sizeof(double)*(1 + params.dimCount);
			char *indCh = new char[sizeI];
			best->indToCh(indCh);
			if (rank == 0){
				t2 = MPI_Wtime();
				Particle temp1, temp2;
				temp1 = *best;
				char *buf = new char[sizeI*size+sizeI];
				MPI_Gather(indCh, sizeI, MPI_CHAR, buf, sizeI*size, MPI_CHAR, 0, MPI_COMM_WORLD);
				for (int i = 0; i < size; ++i){
					temp2.chToInd(buf + i*sizeI);
					if (temp2.fitness < temp1.fitness)temp1 = temp2;
				}
				delete[] buf;
				if (task.check(temp1.fitness, 0.001))++successTests;
				cout << "The best fitness value: " << std::fixed << temp1.fitness << std::endl;
				cout << "Position: [" << temp1.pos[0];
				for (int i = 1; i < params.dimCount; ++i)cout << "," << temp1.pos[i];
				cout << "]" << std::endl;
				cout << "Time: " << t2 - t1 << std::endl;
				cout << "Iteration: " << iteration << std::endl;
			}
			else{
				MPI_Gather(indCh, sizeI, MPI_CHAR, 0, 0, MPI_CHAR, 0, MPI_COMM_WORLD);
			}
			delete[] indCh;
		}
	}
}
