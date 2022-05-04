#ifndef GENETICALGORITHM_KERNEL

#define GENETICALGORITHM_KERNEL

#include <time.h>

typedef struct GENETICALGORITHM {
	int no_Dimension;			//number of dimension
	double **SearchSpace;		//the array to store the search space

	int no_GeneBit;

	int Mutation_Mode;			//different mutation mode
	double Mutation_Rate;

	int CrossOver_Mode;			//different crossover mode
	int no_CrossOver_Point;		//number of crossover point
	double CrossOver_Rate;

	int Selection_Mode;			//different selection mode
	int Selection_TournamentSize;

	int cur_Population_Size;
	int nxt_Population_Size;	// How many chromosomes will be generated

	double **cur_Population;
	double **tmp_Population;
	double **nxt_Population;

	double *cur_Fitness;
	double *tmp_Fitness;
	double *nxt_Fitness;

	double *Optimal_Chromosome;		//pointer pointing to the optimal chromosome
	double Optimal_Fitness;

	int *Chromosome_Index;
	double *Chromosome_Fitness;
	double *Chromosome_SelectionPressure;
	double *Chromosome_Probability;

	double *Convergence;
	double *Diversity;
	double *Processing_Time;

	int no_Generation;				//number of generation
} GeneticAlgorithm;

extern int GENETIC_ALGORITHM_TSP_MUTATION_POINT;
extern int GENETIC_ALGORITHM_TSP_MUTATION_MODE;

#define GA_MAX_GENERATION 50000
#define GA_STOPPING_GENERATION 10

void GeneticAlgorithm_Construction(GeneticAlgorithm*, int, int, int, double, double);
void GeneticAlgorithm_Destruction(GeneticAlgorithm);
void GeneticAlgorithm_Selection_Standard_Elitism(GeneticAlgorithm*);
void GeneticAlgorithm_UnConstraint_CrossOver_Uniform(double*, double*, double**, double**, int, double);
void GeneticAlgorithm_UnConstraint_Mutation_RealCode_Gaussian(double*, double**, double**, int, double, int);
void GeneticAlgorithm_UnConstraint_Standard(GeneticAlgorithm*, int, int, double, double (*fn_ptr)(double*));

#endif
