
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <conio.h>
#include <chrono>
#include <cmath>
const double PI = 3.1415926535897932384626433832795;
using namespace std;
inline double rndDouble(double min, double max){
	return (rand()/ double(RAND_MAX))*(max - min) + min;
}
struct {
	int populationSize;
	int genotypeSize;
	double pMutation;
	double pCrossover;
	int eliteCount;
}params;
struct individual{
	double fitness;
	double *genotype;
	individual(){
		genotype = new double[params.genotypeSize];
	}
	~individual(){
		delete[] genotype;
	}
	individual(const individual& i){
		genotype = new double[params.genotypeSize];
		copy(i.genotype, i.genotype + params.genotypeSize, genotype);
		fitness = i.fitness;
	}
	const individual& operator=(const individual& i){
		copy(i.genotype, i.genotype + params.genotypeSize, genotype);
		fitness = i.fitness;
		return *this;
	}
	void indToCh(char* dest){
		memcpy(dest, &fitness, sizeof(fitness));
		dest += sizeof(fitness);
		memcpy(dest, genotype, sizeof(genotype[0])*params.genotypeSize);
	}
	void chToInd(const char* dest){
		memcpy(&fitness, dest, sizeof(fitness));
		dest += sizeof(fitness);
		memcpy(genotype, dest, sizeof(genotype[0])*params.genotypeSize);
	}
};
double fitness_function(const double *genotype){
	//Rastrigin's function
	/*double temp = 0;
	for (int i = 0; i < params.genotypeSize; ++i)
	temp += (genotype[i] * genotype[i] - 10 * cos(2 * PI*genotype[i]));
	return 10 * params.genotypeSize + temp;*/
	//Schwefel's function
	double temp = 0;
	for (int i = 0; i < params.genotypeSize; ++i)temp += (-genotype[i] * sin(sqrt(abs(genotype[i]))));
	return temp;
	//Rosenbrock's function
	/*double temp = 0;
	for (int i = 0; i < params.genotypeSize - 1; ++i){
		temp += 100 * (genotype[i + 1] - genotype[i] * genotype[i])*(genotype[i + 1] - genotype[i] * genotype[i]) + (genotype[i] - 1)*(genotype[i] - 1);
	}
	return temp;
	double temp1 = 0, temp2 = 1;
	for (int i = 0; i < params.genotypeSize; ++i){
		temp1 += genotype[i] * genotype[i] / 4000.0;
		temp2 *= cos(genotype[i] / sqrt(double(i+1)));
	}
	return 1 + temp1 - temp2;*/
}
double fitness_function(const individual &solution){
	return fitness_function(solution.genotype);
}
bool isEqual(double a, double b, double delta){
	if (abs(a - b) < delta)return true;
	else return false;
}
int main(int argc, char **argv){
	//params initializing
	params.populationSize = 200;
	params.genotypeSize = 20;
	params.pMutation = 0.01;
	params.pCrossover = 0.8;
	params.eliteCount = 0.05*params.populationSize;
	{
		int tempi;
		double tempd;
		bool isElitCountCh = false;
		for (int i = 0; i < argc; ++i){
			if (strncmp(argv[i], "ps=", strlen("ps=")) == 0){
				tempi = atoi(argv[i] + strlen("ps="));
				if (tempi != 0)	params.populationSize = tempi;
			}
			else if (strncmp(argv[i], "gs=", strlen("gs=")) == 0){
				tempi = atoi(argv[i] + strlen("gs="));
				if (tempi != 0) params.genotypeSize = tempi;
			}
			else if (strncmp(argv[i], "cp=", strlen("cp=")) == 0){
				tempd = atof(argv[i] + strlen("cp="));
				if (tempd != 0) params.pCrossover = tempd;
			}
			else if (strncmp(argv[i], "mp=", strlen("mp=")) == 0){
				tempd = atof(argv[i] + strlen("mp="));
				if (tempd != 0) params.pMutation = tempd;
			}
			else if (strncmp(argv[i], "el=", strlen("el=")) == 0){
				tempi = atoi(argv[i] + strlen("el="));
				if (tempi != 0){
					params.eliteCount = tempi;
					isElitCountCh = true;
				}
			}
			if (!isElitCountCh)params.eliteCount = 0.05*params.populationSize;
		}
	}
	srand(time(0));
	double minGenValue = -600, maxGenValue = 600;
	std::vector<individual> population;
	population.reserve(params.populationSize);
	std::vector<individual> parentPool;
	parentPool.reserve(params.populationSize);
	chrono::high_resolution_clock::time_point t1, t2;
	t1 = chrono::high_resolution_clock::now();
	//population generating
	double size = maxGenValue - minGenValue;
	for (int i = 0; i<params.populationSize; ++i){
		individual temp;
		for (int j = 0; j < params.genotypeSize;++j)
			temp.genotype[j] = rndDouble(minGenValue, maxGenValue);
		temp.fitness = fitness_function(temp);
		population.push_back(temp);
	}
	vector<double> individProp(params.populationSize);
	double sum;
	//basic genotype cicle
	int ji = 0;
	double tfw = population[0].fitness;
	int maxIter = 2000,maxRepeat=200,repeatCount=0;
	while (ji<maxIter /*&& repeatCount<maxRepeat*/){
		++ji;
		//forming parent pool
		parentPool.clear();
		
		double maxE = max_element(population.begin(), population.end(), [](individual &a, individual &b){return a.fitness < b.fitness; })->fitness;
		sum = 0;
		std::for_each(population.begin(), population.end(), [&sum](individual &x){sum += -x.fitness; });
		sum += params.populationSize*maxE;
		std::transform(population.begin(), population.end(), individProp.begin(), [&sum,&maxE](individual &x){return (-x.fitness+maxE)/sum; });
		sum = 1;
		//std::for_each(individProp.begin(), individProp.end(), [&sum](double &x){sum += x; });
		for (int i = 0; i<params.populationSize; ++i){
			double randt = rndDouble(0, sum);
			int j = -1;
			double tempSum = 0;
			do{
				++j;
				tempSum += individProp[j];
			} while (tempSum<randt && j<params.populationSize - 1);
			parentPool.push_back(population[j]);
		}
		//forming new generation
		if (params.eliteCount>0){
			sort(population.begin(), population.end(), [](individual& a, individual& b){return a.fitness < b.fitness; });
		}
		population.erase(population.begin() + params.eliteCount, population.end());
		while (population.size()<params.populationSize){
			int firstParent = rand() % params.populationSize;
			int secondParent = rand() % params.populationSize;
			if (rndDouble(0, 1) < params.pCrossover){
				individual child;
				/// BLX -alpha crossover
				/*{
					double alpha = 0.5;
					for (int i = 0; i < genotypeSize; ++i){
						double cmin = min(parentPool[firstParent].genotype[i], parentPool[secondParent].genotype[i]);
						double cmax = max(parentPool[firstParent].genotype[i], parentPool[secondParent].genotype[i]);
						double ii = cmax - cmin;
						child.genotype[i] = rndDouble(cmin - ii*alpha, cmax + ii*alpha);
						if (child.genotype[i] < minGenValue)child.genotype[i] = minGenValue;
						if (child.genotype[i]> maxGenValue)child.genotype[i] = maxGenValue;
					}
					child.fitness = fitness_function(child);
					population.push_back(child);
				}*/
				///SBX crossover
				{
					double betta, u;
					double const n = 2;
					if ((u = rndDouble(0, 1)) > 0.5001)
						betta = pow(double(1) / (2 - 2 * u), double(1) / (n + 1));
					else 
						betta = pow(2 * u, double(1) / (n + 1));
					for (int i = 0; i < params.genotypeSize; ++i){
						child.genotype[i] = 0.5*((1 - betta)*parentPool[firstParent].genotype[i] + (1 + betta)*parentPool[secondParent].genotype[i]);
						if (child.genotype[i] < minGenValue)child.genotype[i] = minGenValue;
						if (child.genotype[i]> maxGenValue)child.genotype[i] = maxGenValue;
					}
					child.fitness = fitness_function(child);
					population.push_back(child);
					if (population.size()<params.populationSize){
						if ((u = rndDouble(0, 1)) > 0.5001)betta = pow(double(1) / (2 - 2 * u), double(1) / (n + 1));
						else betta = pow(2 * u, double(1) / (n + 1));
						for (int i = 0; i < params.genotypeSize; ++i){
							child.genotype[i] = 0.5*((1 - betta)*parentPool[firstParent].genotype[i] + (1 + betta)*parentPool[secondParent].genotype[i]);
							if (child.genotype[i] < minGenValue)child.genotype[i] = minGenValue;
							if (child.genotype[i]> maxGenValue)child.genotype[i] = maxGenValue;
						}
						child.fitness = fitness_function(child);
						population.push_back(child);
					}
				}
				/// Arithmetical crossover
				/*
				double lambda = rand() / double(RAND_MAX);
				for (int i = 0; i < genotypeSize; ++i)
				child.genotype[i] = lambda*parentPool[firstParent].genotype[i] + (1 - lambda)*parentPool[secondParent].genotype[i];
				child.fitness = fitness_function(child);
				population.push_back(child);
				lambda = rand() / double(RAND_MAX);
				for (int i = 0; i < genotypeSize; ++i)
				child.genotype[i] = lambda*parentPool[secondParent].genotype[i] + (1 - lambda)*parentPool[firstParent].genotype[i];
				child.fitness = fitness_function(child);
				population.push_back(child);*/
			}
			else{
				population.push_back(parentPool[firstParent]);
				if (population.size()<params.populationSize)population.push_back(parentPool[secondParent]);
			}
		}
		//uniform mutation
		/*
		for (int i = 0; i < populationSize; ++i){
			if (pMutation>(rand() / double(RAND_MAX)))population[i].genotype[0] = rndDouble(minGenValue, maxGenValue);
			if (pMutation>(rand() / double(RAND_MAX)))population[i].genotype[1] = rndDouble(minGenValue, maxGenValue);
		}*/
		//non-uniform mutation Mihalevich
		double b = 2;
		for (int i = params.eliteCount>0 ? 1 : 0; i < params.populationSize; ++i){
			for (int j = 0; j<params.genotypeSize; ++j){
				if (params.pMutation>(rand() / double(RAND_MAX))){
					if (rndDouble(0, 1)<0.5)
						population[i].genotype[j] += (maxGenValue - population[i].genotype[j])*(1 - pow(rndDouble(0, 1), pow((1 - ji / maxIter), b)));
					else
						population[i].genotype[j]-=(population[i].genotype[j] - minGenValue)*(1 - pow(rndDouble(0, 1), pow((1 - ji / maxIter), b)));
					population[i].fitness = fitness_function(population[i]);
				}
			}
		}
		//output statistic info
		int minI = 0, maxI = 0;
		for (int i = 1; i<params.populationSize; ++i){
			if (population[i].fitness<population[minI].fitness)minI = i;
			if (population[i].fitness>population[maxI].fitness)maxI = i;
		}
		cout <<ji<<":"<< "Min|Max fitness values " << setw(12) << population[minI].fitness <<" | "<< setw(12) << population[maxI].fitness << endl;
		if (isEqual(tfw, population[minI].fitness, 0.001)){
			++repeatCount;
		}
		else{
			repeatCount = 0;
			tfw = population[minI].fitness;
		}
	}
	t2 = chrono::high_resolution_clock::now();
	sort(population.begin(), population.end(), [](individual& a, individual& b){return a.fitness < b.fitness; });
	cout << ji << ":" << "Min|Max fitness values " << setw(12) << population.begin()->fitness << " | " << setw(12) << (population.end()-1)->fitness << endl;
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() / 1000000.0;
	cout << "Time: " << duration<<endl;
	cout << std::fixed << population.begin()->fitness;
	getch();
	return 0;
}


