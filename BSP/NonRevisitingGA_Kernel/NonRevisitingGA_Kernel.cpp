/******************************************************************************************/
//	Program Name	:	Non-Revisiting Genetic Algorithm (NrGA)
//	File Name		:	NonRevisitingGA_Kernel.c
//	Auther			:	Dr. Chow Chi Kin, Dr. Yuen Shiu Yin
//	Edit by			:	Leung Shing Wa
//	University		:	City University of Hong Kong
//	Department		:	Electronic Engineering
//	Last Update		:	10 Sep 2009
//	Reference		:	A Genetic Algorithm that Adaptively Mutates and Never Revisits, 
//						IEEE Transactions on Evolutionary Computation, 
//						Vol 13(2) (April 2009) 454-472.
//
//	Discription		:	The main body of NrGA.
#include "../../BSP/NonRevisitingGA_Kernel/NonRevisitingGA_Kernel.h"

#include<stdio.h>
#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#include<malloc.h>
#include<math.h>
#include<time.h>

#include "../../BSP/ACC_Math_Kernel/ACC_Math_Kernel.h"
#include "../../BSP/RecordRevisitingScheme_Kernel/RecordRevisitingScheme_Kernel.h"
#include "../NonRevisitingScheme_Kernel/NonRevisitingScheme_Kernel.h"
#include "../ClusterTree_Kernel/ClusterTree_Kernel.h"
#include "../GeneticAlgorithm_Kernel/GeneticAlgorithm_Kernel.h"

// build the non-revisiting GA by the giving parameter.
void NonRevisitingGA_Construction(NonRevisiting_GA *cur_GA, int no_dimension, int cur_population_size, int nxt_population_size, double crossover_rate, double mutation_rate, int max_archive_size, int activate)
{
	(*cur_GA).NRS_Info.Activate = activate;
	(*cur_GA).NRS_Info.Max_no_Leaf = 0;

	GeneticAlgorithm_Construction(&(*cur_GA).GA_Info, no_dimension, cur_population_size, nxt_population_size, crossover_rate, mutation_rate);
	Cluster_Tree_Construction(&(*cur_GA).NRS_Info.Solution_Tree, no_dimension, 1, max_archive_size);

	return;
}

// destroy the GA.
void NonRevisitingGA_Destruction(NonRevisiting_GA cur_GA)
{

	GeneticAlgorithm_Destruction(cur_GA.GA_Info);

	Cluster_Tree_Destruction(cur_GA.NRS_Info.Solution_Tree);

	return;
}
/*****************************************
 *  Function: mutate the input chromosome.
 *  		Firstly, it will check whether the chromosome is revisited or not by comparing
 *  		the chromosome with the cluster tree inside the cur_NRS.
 *  		And it will mutate the chromosome if it is revisited.
 *
 *  Input parameters:
 *		cur_Chromosome: the chromosome for mutation
 *		no_gene: the number of gene in the input chromosome
 *		cur_NRS: the input Non revisiting scheme, it provides a cluster tree
 */
void NonRevisitingGA_UnConstraint_Mutation(double **cur_Chromosome, int no_gene, NonRevisitingScheme *cur_NRS)
{
	int i;
	double node_difference;
	Cluster_Node *cur_Node;

	NonRevisitingScheme_UnConstraint_NodeSearch(cur_Chromosome, no_gene, &(*cur_NRS).Solution_Tree, &cur_Node);
	//返回查找到的叶节点，如果访问过则进行变异

	//============ Mutate Gene
	node_difference = Cluster_Node_Difference(cur_Node, *cur_Chromosome, no_gene);
	// 检查查找返回的叶节点与当前节点是否一样，如果一样则node_difference = 0； 否则，没有访问过

	if(node_difference == 0)	//if revisited
	{
		printf("befor mutate:\n");
		for (i = 0; i < (*cur_NRS).Solution_Tree.no_Dimension; ++i)
		{
			printf("[%.1f, %.1f]:%.1f\t", (*cur_NRS).Solution_Tree.cur_Interval[i][0], (*cur_NRS).Solution_Tree.cur_Interval[i][1], (*cur_Chromosome)[i]);
		}
		//printf("\noptimal node:\n");
		//for (i = 0; i < (*cur_NRS).Solution_Tree.no_Dimension; ++i)
		//{
		//	printf("%f\t", (*cur_Node).Optimal_Data[i]);
		//}
		//printf("\n");


		do {
			//逐维检查是否还剩余可以划分的空间
			for(i=0; i < (*cur_NRS).Solution_Tree.no_Dimension; ++i)
			{
				if((*cur_NRS).Solution_Tree.cur_Interval[i][1] - (*cur_NRS).Solution_Tree.cur_Interval[i][0] > 1.0)	 //if there is any space to change
//					(*cur_Chromosome)[i] = (int) ((double) rand() / (RAND_MAX + 1) * ((*cur_NRS).Solution_Tree.cur_Interval[i][1] - (*cur_NRS).Solution_Tree.cur_Interval[i][0])) + (*cur_NRS).Solution_Tree.cur_Interval[i][0] + 0.5;
					//如果有，则变异当前个体。划分当前子空间
					(*cur_Chromosome)[i] = (int) ((double) rand() / (RAND_MAX + 1.0) * ((*cur_NRS).Solution_Tree.cur_Interval[i][1] - (*cur_NRS).Solution_Tree.cur_Interval[i][0])) + (*cur_NRS).Solution_Tree.cur_Interval[i][0] + 0.5;

			}
			// 逐维检查是否有维度变异了。如果变异了则跳出循环
			for (i = 0; i < (*cur_NRS).Solution_Tree.no_Dimension; ++i)
			{
				if ((*cur_Chromosome)[i] != (*cur_Node).Optimal_Data[i])  //Optimal_Data = the evaluated point in this partition
					break;												 //once a difference occur, break the mutation
			}

		} while(i == (*cur_NRS).Solution_Tree.no_Dimension);

		printf("\nafter mutate:\n");
		for (i = 0; i < (*cur_NRS).Solution_Tree.no_Dimension; ++i)
		{
			printf("[%.1f, %.1f]:%.1f\t", (*cur_NRS).Solution_Tree.cur_Interval[i][0], (*cur_NRS).Solution_Tree.cur_Interval[i][1], (*cur_Chromosome)[i]);
		}
		printf("\n");
		system("pause");


	}
		//进行检查，变异产生的值是否超出了边界
		for(i=0; i < (*cur_NRS).Solution_Tree.no_Dimension; ++i)
			if((*cur_Chromosome)[i] < (*cur_NRS).Solution_Tree.Interval[i][0] || (*cur_Chromosome)[i] > (*cur_NRS).Solution_Tree.Interval[i][1])
			{
				printf("%lf %lf %lf Error !!!\n", (*cur_Chromosome)[i], (*cur_NRS).Solution_Tree.cur_Interval[i][0], (*cur_NRS).Solution_Tree.cur_Interval[i][1]);
	//			Pause();
			}

	return;
}
void NonRevisitingGA_UnConstraint_Mutation_V4(double **cur_Chromosome, int no_gene, NonRevisitingScheme *cur_NRS)
{
	int i;
	double node_difference;
	Cluster_Node *cur_Node;


	NonRevisitingScheme_UnConstraint_NodeSearch(cur_Chromosome, no_gene, &(*cur_NRS).Solution_Tree, &cur_Node);
	//返回查找到的叶节点，如果访问过则进行变异

	//============ Mutate Gene
	node_difference = Cluster_Node_Difference(cur_Node, *cur_Chromosome, no_gene);
	// 检查查找返回的叶节点与当前节点是否一样，如果一样则node_difference = 0； 否则，没有访问过

	if (node_difference == 0)	//if revisited
	{
		printf("befor mutate:\n");
		for (i = 0; i < (*cur_NRS).Solution_Tree.no_Dimension; ++i)
		{
			printf("[%.1f, %.1f]:%.1f\t", (*cur_NRS).Solution_Tree.cur_Interval[i][0], (*cur_NRS).Solution_Tree.cur_Interval[i][1], (*cur_Chromosome)[i]);
		}
		//printf("\noptimal node:\n");
		//for (i = 0; i < (*cur_NRS).Solution_Tree.no_Dimension; ++i)
		//{
		//	printf("%f\t", (*cur_Node).Optimal_Data[i]);
		//}
		//printf("\n");


		do {
			//逐维检查是否还剩余可以划分的空间
			for (i = 0; i < (*cur_NRS).Solution_Tree.no_Dimension; ++i)
			{
				if ((*cur_NRS).Solution_Tree.cur_Interval[i][1] - (*cur_NRS).Solution_Tree.cur_Interval[i][0] > 1.0)	 //if there is any space to change
																														 //					(*cur_Chromosome)[i] = (int) ((double) rand() / (RAND_MAX + 1) * ((*cur_NRS).Solution_Tree.cur_Interval[i][1] - (*cur_NRS).Solution_Tree.cur_Interval[i][0])) + (*cur_NRS).Solution_Tree.cur_Interval[i][0] + 0.5;
																														 //如果有，则变异当前个体。划分当前子空间
					(*cur_Chromosome)[i] = (int)((double)rand() / (RAND_MAX + 1.0) * ((*cur_NRS).Solution_Tree.cur_Interval[i][1] - (*cur_NRS).Solution_Tree.cur_Interval[i][0])) + (*cur_NRS).Solution_Tree.cur_Interval[i][0] + 0.5;

			}
			// 逐维检查是否有维度变异了。如果变异了则跳出循环
			for (i = 0; i < (*cur_NRS).Solution_Tree.no_Dimension; ++i)
			{
				if ((*cur_Chromosome)[i] != (*cur_Node).Optimal_Data[i])  //Optimal_Data = the evaluated point in this partition
					break;												 //once a difference occur, break the mutation
			}

		} while (i == (*cur_NRS).Solution_Tree.no_Dimension);

		printf("\nafter mutate:\n");
		for (i = 0; i < (*cur_NRS).Solution_Tree.no_Dimension; ++i)
		{
			printf("[%.1f, %.1f]:%.1f\t", (*cur_NRS).Solution_Tree.cur_Interval[i][0], (*cur_NRS).Solution_Tree.cur_Interval[i][1], (*cur_Chromosome)[i]);
		}
		printf("\n");
		system("pause");


	}
	//进行检查，变异产生的值是否超出了边界
	for (i = 0; i < (*cur_NRS).Solution_Tree.no_Dimension; ++i)
		if ((*cur_Chromosome)[i] < (*cur_NRS).Solution_Tree.Interval[i][0] || (*cur_Chromosome)[i] > (*cur_NRS).Solution_Tree.Interval[i][1])
		{
			printf("%lf %lf %lf Error !!!\n", (*cur_Chromosome)[i], (*cur_NRS).Solution_Tree.cur_Interval[i][0], (*cur_NRS).Solution_Tree.cur_Interval[i][1]);
			//			Pause();
		}

	return;
}

void NonRevisitingGA_UnConstraint_Standard(NonRevisiting_GA *cur_GA, double optimization_mode, int terminate_mode, double terminate_parameter, int tree_mode, double (*fn_ptr)(double *chromosome))
{
	//============ terminate_mode: 0 = Fixed Generation    (terminate_parameter = max. no. of iteration)
	//                             1 = Fixed Accuracy	   (terminate_parameter = target accuracy)
	//                             2 = Fitness Improvement (terminate_parameter = delta imrpovement)
	//
	//============ optimization_mode: 0 = Minization, 1: Maximization

	int i, j;
	int chromosome_index_1, chromosome_index_2, gene_index;
	int optimal_chromosome_index;
	long clock_begin;

	int no_stopping_generation;
	double prev_optimum;

	char GA_filename[256];
	char GA_filename2[256];
	int cluster_write_mode;
	FILE *GA_ptr = NULL;

	cluster_write_mode = 1;
	if(cluster_write_mode == 1)
	{
		GA_ptr = fopen("GA_Solution", "w");
		fprintf(GA_ptr, "%d %d\n", (int) terminate_parameter, (*cur_GA).GA_Info.cur_Population_Size);
	}

	clock_begin = clock();

	no_stopping_generation = 0;
	(*cur_GA).GA_Info.no_Generation = 0;

	//============ Initialize Chromosome Tree
	if(tree_mode == 0)
	{
		(*cur_GA).NRS_Info.Max_no_Leaf = 0;
		for(i=0; i< (*cur_GA).GA_Info.no_Dimension; ++i)
			for(j=0; j< 2; ++j)
				(*cur_GA).NRS_Info.Solution_Tree.Interval[i][j] = (*cur_GA).GA_Info.SearchSpace[i][j];
	}

	//============ Population Initialization
	for(i=0;i< (*cur_GA).GA_Info.cur_Population_Size; ++i)
	{
		for(j=0;j< (*cur_GA).GA_Info.no_Dimension; ++j)
			(*cur_GA).GA_Info.cur_Population[i][j] = rand() % (int) ((*cur_GA).GA_Info.SearchSpace[j][1] - (*cur_GA).GA_Info.SearchSpace[j][0]) + (*cur_GA).GA_Info.SearchSpace[j][0] + 0.5;

		if((*cur_GA).NRS_Info.Activate == 1 && (*(*cur_GA).NRS_Info.Solution_Tree.Root).Dimension != -2)
			NonRevisitingGA_UnConstraint_Mutation(&(*cur_GA).GA_Info.cur_Population[i], (*cur_GA).GA_Info.no_Dimension, &(*cur_GA).NRS_Info);
		else
			printf("root dimension = -2/n");
		//评估适应值，默认为最大化问题。如果是最小化问题需要添加-1
		(*cur_GA).GA_Info.cur_Fitness[i] = (1.0 - 2.0 * optimization_mode) * (*fn_ptr)((*cur_GA).GA_Info.cur_Population[i]);

		//V2 记录点
		//RecordRevisiting((*cur_GA).GA_Info.cur_Population[i], (*cur_GA).GA_Info.cur_Fitness[i], (*cur_GA).GA_Info.no_Dimension,
		//		&(*cur_GA).NRS_Info, &(*cur_GA).NRS_Info.Solution_Tree.Current_ID, (*cur_GA).NRS_Info.Solution_Tree.cur_Interval);
		NonRevisitingScheme_UnConstraint_Update(&(*cur_GA).NRS_Info, (*cur_GA).GA_Info.cur_Population[i], (*cur_GA).GA_Info.cur_Fitness[i]);

		// 对于只寻找全局最优的情况，直接更新记录
		if(i==0 || (*cur_GA).GA_Info.cur_Fitness[i] < (*cur_GA).GA_Info.Optimal_Fitness)
		{
			(*cur_GA).GA_Info.Optimal_Fitness = (*cur_GA).GA_Info.cur_Fitness[i];
			optimal_chromosome_index = i;
		}

		//===================== BSP Tree Investigation =================================
		if(cluster_write_mode == 1)
		{
			for(j=0;j< (*cur_GA).GA_Info.no_Dimension; ++j)
				fprintf(GA_ptr, "%lf ", (*cur_GA).GA_Info.cur_Population[i][j]);
			fprintf(GA_ptr, "\n");
		}
		//==============================================================================
	}

	for(j=0;j< (*cur_GA).GA_Info.no_Dimension; ++j)
		(*cur_GA).GA_Info.Optimal_Chromosome[j] = (*cur_GA).GA_Info.cur_Population[optimal_chromosome_index][j];

	++(*cur_GA).GA_Info.no_Generation;

	//===================== BSP Tree Investigation =================================
	if(cluster_write_mode == 1)
	{
		sprintf(GA_filename, "GA_Tree_%d", (*cur_GA).GA_Info.no_Generation);
//		printf("%s\n", GA_filename);
		Cluster_Tree_Write(&(*cur_GA).NRS_Info.Solution_Tree, GA_filename);

		sprintf(GA_filename2, "GA_Tree_%d_noPrune", (*cur_GA).GA_Info.no_Generation);
//		printf("%s\n", GA_filename2);
		Record_Cluster_Tree_Write(&(*cur_GA).NRS_Info.Solution_Tree, GA_filename2);
	}
	//==============================================================================

	(*cur_GA).GA_Info.Convergence = (double *)realloc((*cur_GA).GA_Info.Convergence, sizeof(double) * (*cur_GA).GA_Info.no_Generation);
	(*cur_GA).GA_Info.Convergence[(*cur_GA).GA_Info.no_Generation - 1] = (*cur_GA).GA_Info.Optimal_Fitness;

	(*cur_GA).GA_Info.Processing_Time = (double *)realloc((*cur_GA).GA_Info.Processing_Time, sizeof(double) * (*cur_GA).GA_Info.no_Generation);
	(*cur_GA).GA_Info.Processing_Time[(*cur_GA).GA_Info.no_Generation - 1] = (double) (clock() - clock_begin) / CLOCKS_PER_SEC;
	//===================================================================================================

	//============ Evoluation
	do {
		//============ Reproduction
		//------------ Cross-Over

		for(i=0; i< (*cur_GA).GA_Info.nxt_Population_Size / 2; ++i)
		{
			chromosome_index_1 = rand() % (*cur_GA).GA_Info.cur_Population_Size;
			chromosome_index_2 = rand() % (*cur_GA).GA_Info.cur_Population_Size;
			while(chromosome_index_2 == chromosome_index_1)
				chromosome_index_2 = rand() % (*cur_GA).GA_Info.cur_Population_Size;


			GeneticAlgorithm_UnConstraint_CrossOver_Uniform((*cur_GA).GA_Info.cur_Population[chromosome_index_1],
																	(*cur_GA).GA_Info.cur_Population[chromosome_index_2],
																	&(*cur_GA).GA_Info.nxt_Population[i*2],
																	&(*cur_GA).GA_Info.nxt_Population[i*2+1],
																	(*cur_GA).GA_Info.no_Dimension,
																	(*cur_GA).GA_Info.CrossOver_Rate);
		}

		//------------ Mutation
		for(i=0; i< (*cur_GA).GA_Info.nxt_Population_Size; ++i)
		{

			if((*cur_GA).NRS_Info.Activate == 0)
				GeneticAlgorithm_UnConstraint_Mutation_RealCode_Gaussian((*cur_GA).GA_Info.nxt_Population[i],
																		 &(*cur_GA).GA_Info.nxt_Population[i],
																		 (*cur_GA).GA_Info.SearchSpace,
																		 (*cur_GA).GA_Info.no_Dimension,
																		 (*cur_GA).GA_Info.Mutation_Rate, 0);

			else
			{
				gene_index = rand() % (*cur_GA).GA_Info.no_Dimension;

				(*cur_GA).GA_Info.nxt_Population[i][gene_index] += (2 * (rand() % 2) - 1);			// to generate +1/-1 disturbing
				if((*cur_GA).GA_Info.nxt_Population[i][gene_index] < (*cur_GA).GA_Info.SearchSpace[gene_index][0] + 0.5)
					(*cur_GA).GA_Info.nxt_Population[i][gene_index] = (*cur_GA).GA_Info.SearchSpace[gene_index][0] + 0.5;
				else if((*cur_GA).GA_Info.nxt_Population[i][gene_index] > (*cur_GA).GA_Info.SearchSpace[gene_index][1] - 0.5)
					(*cur_GA).GA_Info.nxt_Population[i][gene_index] = (*cur_GA).GA_Info.SearchSpace[gene_index][1] - 0.5;

			}

			if((*cur_GA).NRS_Info.Activate == 1 && (*(*cur_GA).NRS_Info.Solution_Tree.Root).Dimension != -2)
				NonRevisitingGA_UnConstraint_Mutation(&(*cur_GA).GA_Info.nxt_Population[i], (*cur_GA).GA_Info.no_Dimension, &(*cur_GA).NRS_Info);

			(*cur_GA).GA_Info.nxt_Fitness[i] = (1.0 - 2.0 * optimization_mode) * (*fn_ptr)((*cur_GA).GA_Info.nxt_Population[i]);

			NonRevisitingScheme_UnConstraint_Update(&(*cur_GA).NRS_Info, (*cur_GA).GA_Info.nxt_Population[i], (*cur_GA).GA_Info.nxt_Fitness[i]);
			//RecordRevisiting((*cur_GA).GA_Info.nxt_Population[i], (*cur_GA).GA_Info.nxt_Fitness[i], (*cur_GA).GA_Info.no_Dimension,
			//		&(*cur_GA).NRS_Info, &(*cur_GA).NRS_Info.Solution_Tree.Current_ID, (*cur_GA).NRS_Info.Solution_Tree.cur_Interval);


			//===================== BSP Tree Investigation =================================
			if(cluster_write_mode == 1)
			{
				for(j=0;j< (*cur_GA).GA_Info.no_Dimension; ++j)
					fprintf(GA_ptr, "%lf ", (*cur_GA).GA_Info.nxt_Population[i][j]);
				fprintf(GA_ptr, "\n");
			}
			//==============================================================================
		}

	
		//------------ Selection
		GeneticAlgorithm_Selection_Standard_Elitism(&(*cur_GA).GA_Info);

		if(terminate_mode == 2)
		{
			if((*cur_GA).GA_Info.no_Generation > 0)
			{
				if((prev_optimum - (*cur_GA).GA_Info.Optimal_Fitness) / prev_optimum < terminate_parameter)
					++no_stopping_generation;
				else
					no_stopping_generation = 0;
			}
			prev_optimum = (*cur_GA).GA_Info.Optimal_Fitness;
		}

		++(*cur_GA).GA_Info.no_Generation;

		(*cur_GA).GA_Info.Convergence = (double *)realloc((*cur_GA).GA_Info.Convergence, sizeof(double) * (*cur_GA).GA_Info.no_Generation);
		(*cur_GA).GA_Info.Convergence[(*cur_GA).GA_Info.no_Generation - 1] = (*cur_GA).GA_Info.Optimal_Fitness;

		(*cur_GA).GA_Info.Processing_Time = (double *)realloc((*cur_GA).GA_Info.Processing_Time, sizeof(double) * (*cur_GA).GA_Info.no_Generation);
		(*cur_GA).GA_Info.Processing_Time[(*cur_GA).GA_Info.no_Generation - 1] = (double) (clock() - clock_begin) / CLOCKS_PER_SEC;


		if((*cur_GA).GA_Info.no_Generation % 10 == 0)
		{
			printf("Iteration %d: %lf | %d %d\n", (*cur_GA).GA_Info.no_Generation,
												  (*cur_GA).GA_Info.Optimal_Fitness,
												  (*cur_GA).NRS_Info.Solution_Tree.no_Leaf,
												  (*cur_GA).NRS_Info.Max_no_Leaf);
		}

		//===================== BSP Tree Investigation =================================
		if(cluster_write_mode == 1)
		{
			sprintf(GA_filename, "GA_Tree_%d", (*cur_GA).GA_Info.no_Generation);
//			printf("%s\n", GA_filename);
			Cluster_Tree_Write(&(*cur_GA).NRS_Info.Solution_Tree, GA_filename);

			sprintf(GA_filename2, "GA_Tree_%d_noPrune", (*cur_GA).GA_Info.no_Generation);
//			printf("%s\n", GA_filename2);
			Record_Cluster_Tree_Write(&(*cur_GA).NRS_Info.Solution_Tree, GA_filename2);
		}
		//==============================================================================
	} while((terminate_mode == 0 && terminate_parameter > (*cur_GA).GA_Info.no_Generation) ||
			(terminate_mode == 1 && ((1.0 - 2.0 * optimization_mode) * terminate_parameter < (*cur_GA).GA_Info.Optimal_Fitness) && (*cur_GA).GA_Info.no_Generation < NONREVISITINGGA_MAX_GENERATION) ||
			(terminate_mode == 2 && no_stopping_generation < NONREVISITINGGA_STOPPING_GENERATION));



	optimal_chromosome_index = -1;
	for(i=0; i< (*cur_GA).GA_Info.cur_Population_Size; ++i)
	{
		(*cur_GA).GA_Info.cur_Fitness[i] *= (1.0 - 2.0 * optimization_mode);
		if(optimal_chromosome_index == -1 || (*cur_GA).GA_Info.Optimal_Fitness >= (*cur_GA).GA_Info.cur_Fitness[i])
		{
			optimal_chromosome_index = i;
			(*cur_GA).GA_Info.Optimal_Fitness = (*cur_GA).GA_Info.cur_Fitness[i];
		}
	}

	printf("The best Chromosome in this trial is:\n");
	(*cur_GA).GA_Info.Optimal_Fitness *= (1.0 - 2.0 * optimization_mode);
	for(i=0; i< (*cur_GA).GA_Info.no_Dimension; ++i)
	{
		(*cur_GA).GA_Info.Optimal_Chromosome[i] = (*cur_GA).GA_Info.cur_Population[optimal_chromosome_index][i];
		printf("%f ",(*cur_GA).GA_Info.Optimal_Chromosome[i]);
	}
	printf("\n");

	if(cluster_write_mode == 1)
		fclose(GA_ptr);

	return;
}

