/******************************************************************************************/
//	Program Name	:	Non-Revisiting Genetic Algorithm (NrGA)
//	File Name		:	GeneticAlgorithm_Kernel.c
//	Auther			:	Dr. Chow Chi Kin, Dr. Yuen Shiu Yin
//	Edit by			:	Leung Shing Wa
//	University		:	City University of Hong Kong
//	Department		:	Electronic Engineering
//	Last Update		:	10 Sep 2009
//	Reference		:	A Genetic Algorithm that Adaptively Mutates and Never Revisits, 
//						IEEE Transactions on Evolutionary Computation, 
//						Vol 13(2) (April 2009) 454-472.
//
//	Discription		:	Body of the GA.
#include "../../BSP/GeneticAlgorithm_Kernel/GeneticAlgorithm_Kernel.h"

#include<stdio.h>
#include<stdlib.h>
#include<malloc.h>
#include<math.h>
#include<time.h>

#include "../../BSP/ACC_Math_Kernel/ACC_Math_Kernel.h"

int GeneticAlgorithm_Mutation_Mask[256];
int GENETIC_ALGORITHM_TSP_MUTATION_POINT;
int GENETIC_ALGORITHM_TSP_MUTATION_MODE;
//=========function to build the genetic algorithm module
void GeneticAlgorithm_Construction(GeneticAlgorithm *cur_GA, int no_dimension, int cur_population_size, int nxt_population_size, double crossover_rate, double mutation_rate)
{
	int i;

	(*cur_GA).no_Dimension = no_dimension;
	//=========build the search space, each dimension have range
	(*cur_GA).SearchSpace = (double **)malloc(sizeof(double*) * no_dimension);
	for(i=0; i<no_dimension; ++i)
		(*cur_GA).SearchSpace[i] = (double *)malloc(sizeof(double) * 2);  // 2*double, 1 double for max, 1 double for min

	(*cur_GA).Mutation_Rate = mutation_rate;
	(*cur_GA).CrossOver_Rate = crossover_rate;

	(*cur_GA).no_CrossOver_Point = 0;

	(*cur_GA).cur_Population_Size = cur_population_size;
	(*cur_GA).nxt_Population_Size = nxt_population_size;

	(*cur_GA).cur_Population = (double **)malloc(sizeof(double*) * cur_population_size);
	for(i=0; i< (*cur_GA).cur_Population_Size; ++i)
		(*cur_GA).cur_Population[i] = (double *)malloc(sizeof(double) * no_dimension);

	(*cur_GA).tmp_Population = (double **)malloc(sizeof(double*) * cur_population_size);
	for(i=0; i< (*cur_GA).cur_Population_Size; ++i)
		(*cur_GA).tmp_Population[i] = (double *)malloc(sizeof(double) * no_dimension);

	(*cur_GA).nxt_Population = (double **)malloc(sizeof(double*) * nxt_population_size);
	for(i=0; i< (*cur_GA).nxt_Population_Size; ++i)
		(*cur_GA).nxt_Population[i] = (double *)malloc(sizeof(double) * no_dimension);

	(*cur_GA).cur_Fitness = (double *)malloc(sizeof(double) * cur_population_size);
	(*cur_GA).tmp_Fitness = (double *)malloc(sizeof(double) * cur_population_size);
	(*cur_GA).nxt_Fitness = (double *)malloc(sizeof(double) * nxt_population_size);

	(*cur_GA).Optimal_Chromosome = (double *)malloc(sizeof(double) * no_dimension);

	(*cur_GA).Chromosome_Index = (int *)malloc(sizeof(int) * (cur_population_size + nxt_population_size));
	(*cur_GA).Chromosome_Fitness = (double *)malloc(sizeof(double) * (cur_population_size + nxt_population_size));
	(*cur_GA).Chromosome_SelectionPressure = (double *)malloc(sizeof(double) * (cur_population_size + nxt_population_size));
	(*cur_GA).Chromosome_Probability = (double *)malloc(sizeof(double) * (cur_population_size + nxt_population_size));

	(*cur_GA).Processing_Time = (double *)malloc(sizeof(double) * GA_MAX_GENERATION);
	(*cur_GA).Convergence = (double *)malloc(sizeof(double) * GA_MAX_GENERATION);
	(*cur_GA).Diversity = (double *)malloc(sizeof(double) * GA_MAX_GENERATION);

	return;
}



void GeneticAlgorithm_Destruction(GeneticAlgorithm cur_GA)
{
	int i;

	for(i=0; i<cur_GA.no_Dimension; ++i)
		free(cur_GA.SearchSpace[i]);
	free(cur_GA.SearchSpace);

	for(i=0; i< cur_GA.cur_Population_Size; ++i)
		free(cur_GA.cur_Population[i]);
	free(cur_GA.cur_Population);

	for(i=0; i< cur_GA.cur_Population_Size; ++i)
		free(cur_GA.tmp_Population[i]);
	free(cur_GA.tmp_Population);

	for(i=0; i< cur_GA.nxt_Population_Size; ++i)
		free(cur_GA.nxt_Population[i]);
	free(cur_GA.nxt_Population);

	free(cur_GA.cur_Fitness);
	free(cur_GA.tmp_Fitness);
	free(cur_GA.nxt_Fitness);

	free(cur_GA.Optimal_Chromosome);

	free(cur_GA.Chromosome_Index);
	free(cur_GA.Chromosome_Fitness);
	free(cur_GA.Chromosome_SelectionPressure);
	free(cur_GA.Chromosome_Probability);

	free(cur_GA.Processing_Time);
	free(cur_GA.Convergence);
	free(cur_GA.Diversity);

	return;
}


//==========function to uniform unconstraint crossover, each dimension consider once
void GeneticAlgorithm_UnConstraint_CrossOver_Uniform(double *Parent_1, double *Parent_2, double **Child_1, double **Child_2, int no_Dimension, double crossover_rate)
{
	int i;

	for(i=0; i< no_Dimension; ++i)
		if((double) rand() / RAND_MAX <= crossover_rate)
		{
			(*Child_1)[i] = Parent_1[i];
			(*Child_2)[i] = Parent_2[i];
		}
		else
		{
			(*Child_1)[i] = Parent_2[i];
			(*Child_2)[i] = Parent_1[i];
		}

	return;
}

//===========functions to Unconstraint mutation
//mode: 0 = integer, 1 = double
//Mutation rate = variance step size
void GeneticAlgorithm_UnConstraint_Mutation_RealCode_Uniform(double *Parent, double **Child, double **SearchSpace, int no_Dimension, double Mutation_Rate, int mode)
{
	int i;
	double cur_gene;

	for(i=0; i< no_Dimension; ++i)
	{
		cur_gene = Parent[i];
		(*Child)[i] = cur_gene + 2.0 * ((double) (rand() % 2) - 0.5) * ((double) rand() / RAND_MAX * Mutation_Rate) * (SearchSpace[i][1] - SearchSpace[i][0]);
		if(mode == 0)
			(*Child)[i] = (int) (*Child)[i];
	}

	return;
}
// Gaussian distribution to mutate the parents to other place in the search space with the variance step size.
void GeneticAlgorithm_UnConstraint_Mutation_RealCode_Gaussian(double *Parent, double **Child, double **SearchSpace, int no_Dimension, double Mutation_Rate, int mode)
{
	int i;
	double cur_gene;

	for(i=0; i< no_Dimension; ++i)
	{
		cur_gene = Parent[i];
		(*Child)[i] = cur_gene + randn(0.0, Mutation_Rate) * (SearchSpace[i][1] - SearchSpace[i][0]);
		if(mode == 0)
			(*Child)[i] = (int) (*Child)[i];
	}

	return;
}



//======== Standard Elitism selection
void GeneticAlgorithm_Selection_Standard_Elitism(GeneticAlgorithm *GA_Info)
{
	int i, j;

	//======== copy the fitness and index of parents into the array
	for(i=0; i< (*GA_Info).cur_Population_Size; ++i)
	{
		(*GA_Info).Chromosome_Fitness[i] = (*GA_Info).cur_Fitness[i];
		(*GA_Info).Chromosome_Index[i] = i;
	}
	//======== copy the fitness and index of childs into the array
	for(i=0; i< (*GA_Info).nxt_Population_Size; ++i)
	{
		(*GA_Info).Chromosome_Fitness[i + (*GA_Info).cur_Population_Size] = (*GA_Info).nxt_Fitness[i];
		(*GA_Info).Chromosome_Index[i + (*GA_Info).cur_Population_Size] = i + (*GA_Info).cur_Population_Size;
	}

	//======== sorting by fitness
	Sorting_Indexed(&(*GA_Info).Chromosome_Fitness, &(*GA_Info).Chromosome_Index, (*GA_Info).cur_Population_Size + (*GA_Info).nxt_Population_Size);
	for(i=0; i< (*GA_Info).cur_Population_Size; ++i)
	{
		if((*GA_Info).Chromosome_Index[i] < (*GA_Info).cur_Population_Size)
		{
			for(j=0; j < (*GA_Info).no_Dimension; ++j)
				(*GA_Info).tmp_Population[i][j] = (*GA_Info).cur_Population[(*GA_Info).Chromosome_Index[i]][j];
			(*GA_Info).tmp_Fitness[i] = (*GA_Info).Chromosome_Fitness[i];
		}
		else
		{
			for(j=0; j < (*GA_Info).no_Dimension; ++j)
				(*GA_Info).tmp_Population[i][j] = (*GA_Info).nxt_Population[(*GA_Info).Chromosome_Index[i] - (*GA_Info).cur_Population_Size][j];
			(*GA_Info).tmp_Fitness[i] = (*GA_Info).Chromosome_Fitness[i];
		}
	}
	//========= only get the chomosome ranking within the population size
	for(i=0; i< (*GA_Info).cur_Population_Size; ++i)
	{
		for(j=0; j < (*GA_Info).no_Dimension; ++j)
			(*GA_Info).cur_Population[i][j] = (*GA_Info).tmp_Population[i][j];
		(*GA_Info).cur_Fitness[i] = (*GA_Info).tmp_Fitness[i];
		
		//=======find the chomosome with highest fitness
		if((*GA_Info).Optimal_Fitness > (*GA_Info).tmp_Fitness[i])
		{
			for(j=0; j < (*GA_Info).no_Dimension; ++j)
				(*GA_Info).Optimal_Chromosome[j] = (*GA_Info).cur_Population[i][j];
			(*GA_Info).Optimal_Fitness = (*GA_Info).cur_Fitness[i];
		}
	}

	return;
}

// Standard Genetic Algorithm.
void GeneticAlgorithm_UnConstraint_Standard(GeneticAlgorithm *cur_GA, int optimization_mode, int terminate_mode, double terminate_parameter, double (*fn_ptr)(double *chromosome))
{
	//============ terminate_mode: 0 = Fixed Generation    (terminate_parameter = max. no. of iteration)
	//                             1 = Fixed Accuracy	   (terminate_parameter = target accuracy)
	//                             2 = Fitness Improvement (terminate_parameter = delta imrpovement)
	//
	//============ optimization_mode: 0 = Minization, 1: Maximization

	int i, j;
	int chromosome_index_1, chromosome_index_2;
	int optimal_chromosome_index;
	long clock_begin;

	int no_stopping_generation;
	double prev_optimum;

	no_stopping_generation = 0;
	(*cur_GA).no_Generation = 0;
	clock_begin = clock();

	//============ Population Initialization
	for(i=0;i< (*cur_GA).cur_Population_Size; ++i)
	{
		for(j=0;j< (*cur_GA).no_Dimension; ++j)
			(*cur_GA).cur_Population[i][j] = (double) rand() / RAND_MAX * ((*cur_GA).SearchSpace[j][1] - (*cur_GA).SearchSpace[j][0]) + (*cur_GA).SearchSpace[j][0] + 0.5;
		(*cur_GA).cur_Fitness[i] = (1.0 - 2.0 * optimization_mode) * (*fn_ptr)((*cur_GA).cur_Population[i]);

		if(i==0 || (*cur_GA).cur_Fitness[i] < (*cur_GA).Optimal_Fitness)
		{
			(*cur_GA).Optimal_Fitness = (*cur_GA).cur_Fitness[i];
			optimal_chromosome_index = i;
		}
	}

	for(j=0;j< (*cur_GA).no_Dimension; ++j)
		(*cur_GA).Optimal_Chromosome[j] = (*cur_GA).cur_Population[optimal_chromosome_index][j];

	++(*cur_GA).no_Generation;

	(*cur_GA).Convergence = (double *)realloc((*cur_GA).Convergence, sizeof(double) * (*cur_GA).no_Generation);
	(*cur_GA).Convergence[(*cur_GA).no_Generation - 1] = (*cur_GA).Optimal_Fitness;

	(*cur_GA).Processing_Time = (double *)realloc((*cur_GA).Processing_Time, sizeof(double) * (*cur_GA).no_Generation);
	(*cur_GA).Processing_Time[(*cur_GA).no_Generation - 1] = (double) (clock() - clock_begin) / CLOCKS_PER_SEC;
	//===================================================================================================

	//============ Evoluation
	do {
		//============ Reproduction
		//------------ Far Search
		for(i=0; i< (*cur_GA).nxt_Population_Size / 2; ++i)
		{
			chromosome_index_1 = rand() % (*cur_GA).cur_Population_Size;
			chromosome_index_2 = rand() % (*cur_GA).cur_Population_Size;
			while(chromosome_index_2 == chromosome_index_1)
				chromosome_index_2 = rand() % (*cur_GA).cur_Population_Size;


			GeneticAlgorithm_UnConstraint_CrossOver_Uniform((*cur_GA).cur_Population[chromosome_index_1],
															(*cur_GA).cur_Population[chromosome_index_2],
															&(*cur_GA).nxt_Population[i*2],
															&(*cur_GA).nxt_Population[i*2+1],
															(*cur_GA).no_Dimension,
															(*cur_GA).CrossOver_Rate);

		}

		//------------ Local Search
		for(i=0; i< (*cur_GA).nxt_Population_Size; ++i)
		{
/*			GeneticAlgorithm_UnConstraint_Mutation_RealCode_Gaussian((*cur_GA).nxt_Population[i],
																	 &(*cur_GA).nxt_Population[i],
																	 (*cur_GA).SearchSpace,
																	 (*cur_GA).no_Dimension,
																	 (*cur_GA).Mutation_Rate, 1);*/
			GeneticAlgorithm_UnConstraint_Mutation_RealCode_Uniform((*cur_GA).nxt_Population[i],
																	&(*cur_GA).nxt_Population[i],
																	(*cur_GA).SearchSpace,
																	(*cur_GA).no_Dimension,
																	(*cur_GA).Mutation_Rate, 1);
		}

		//------------ Fitness Measurement
		for(i=0; i< (*cur_GA).nxt_Population_Size; ++i)
			(*cur_GA).nxt_Fitness[i] = (1.0 - 2.0 * optimization_mode) * (*fn_ptr)((*cur_GA).nxt_Population[i]);

		//------------ Selection
		GeneticAlgorithm_Selection_Standard_Elitism(&(*cur_GA));

		if(terminate_mode == 2)
		{
			if((*cur_GA).no_Generation > 0)
			{
				if((prev_optimum - (*cur_GA).Optimal_Fitness) / prev_optimum < terminate_parameter)
					++no_stopping_generation;
				else
					no_stopping_generation = 0;
			}
			prev_optimum = (*cur_GA).Optimal_Fitness;
		}

		++(*cur_GA).no_Generation;

		(*cur_GA).Convergence = (double *)realloc((*cur_GA).Convergence, sizeof(double) * (*cur_GA).no_Generation);
		(*cur_GA).Convergence[(*cur_GA).no_Generation - 1] = (*cur_GA).Optimal_Fitness;

		(*cur_GA).Processing_Time = (double *)realloc((*cur_GA).Processing_Time, sizeof(double) * (*cur_GA).no_Generation);
		(*cur_GA).Processing_Time[(*cur_GA).no_Generation - 1] = (double) (clock() - clock_begin) / CLOCKS_PER_SEC;

	} while((terminate_mode == 0 && terminate_parameter > (*cur_GA).no_Generation) ||
			(terminate_mode == 1 && ((1.0 - 2.0 * optimization_mode) * terminate_parameter < (*cur_GA).Optimal_Fitness) && (*cur_GA).no_Generation < GA_MAX_GENERATION) ||
			(terminate_mode == 2 && no_stopping_generation < GA_STOPPING_GENERATION));

	optimal_chromosome_index = -1;
	for(i=0; i< (*cur_GA).cur_Population_Size; ++i)
	{
		(*cur_GA).cur_Fitness[i] *= (1.0 - 2.0 * optimization_mode);
		if(optimal_chromosome_index == -1 || (*cur_GA).Optimal_Fitness > (*cur_GA).cur_Fitness[i])
		{
			optimal_chromosome_index = i;
			(*cur_GA).Optimal_Fitness = (*cur_GA).cur_Fitness[i];
		}
	}

	(*cur_GA).Optimal_Fitness *= (1.0 - 2.0 * optimization_mode);
	for(i=0; i< (*cur_GA).no_Dimension; ++i)
		(*cur_GA).Optimal_Chromosome[i] = (*cur_GA).cur_Population[optimal_chromosome_index][i];

	return;
}

