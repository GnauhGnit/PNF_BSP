/*
 * VisitStrategy.cpp
 *
 *  Created on: Aug 1, 2017
 *      Author: hting
 */


#include "./VisitStrategy.h"

#include <vector>
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/cauchy_distribution.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real.hpp>


int Mutate_Probability;
int Mutate_On_Axis_Range;
int Mutate_To_Nearest_Space;
int Max_No_Visit;
int Converge_Times;
int Mutate_BSP_FULL;
int VisitStrateg_Record(double *cur_Chromosome, double cur_Chromonsome_Fitness, int index,
	NonRevisiting_GA *cur_GA, int *Current_ID, double **cur_Interval, vector<Selection_Strategy> &Select_Vector)
{
	RecordRevisiting((*cur_GA).GA_Info.cur_Population[index], cur_Chromonsome_Fitness, (*cur_GA).GA_Info.no_Dimension,
		&(*cur_GA).NRS_Info, Current_ID, cur_Interval, Select_Vector);
	return 0;
}
void VisitStrategy_On_AxisRange(double *cur_Chromosome, double *tmp_Chromosome, int dim,
	boost::mt19937 &generator, int axis_range)
{

	boost::uniform_real<> uniform_real_generate_r1(0.0, 0.5);
	boost::variate_generator< boost::mt19937&, boost::uniform_real<> > random_real_num_r1(generator, uniform_real_generate_r1);
	boost::uniform_real<> uniform_real_generate_r2(-0.5, 0.0);
	boost::variate_generator< boost::mt19937&, boost::uniform_real<> > random_real_num_r2(generator, uniform_real_generate_r2);
	boost::uniform_real<> uniform_real_generate_r3(-0.5, 0.5);
	boost::variate_generator< boost::mt19937&, boost::uniform_real<> > random_real_num_r3(generator, uniform_real_generate_r3);
	double tmp_r = 0.0;
	for (int j = 0; j < dim; j++)
	{

		if (cur_Chromosome[j] == 0)
		{// 最左端
			tmp_Chromosome[j] = random_real_num_r1() + cur_Chromosome[j];

		}
		else if (cur_Chromosome[j] == axis_range - 1)
		{// 最右端
			tmp_Chromosome[j] = random_real_num_r2() + cur_Chromosome[j];
		}
		else
		{// 中间
			tmp_Chromosome[j] = random_real_num_r3() + cur_Chromosome[j];
		}
		//cout <<  "Mutate on axis: " << cur_Chromosome[j] << "\t|\t" << tmp_Chromosome[j]<<"\tr:"<<tmp_r << endl;
	}
}	// void VisitStrategy_On_AxisRange

void NonRevisitingGA_UnConstraint_Mutation_V4(double **cur_Chromosome, int no_gene, NonRevisitingScheme *cur_NRS, boost::mt19937 &generator)
{
	int i;
	double min_bound, max_bound;
	Cluster_Node *cur_Node = (*cur_NRS).Solution_Tree.Root;
	boost::uniform_real<> uniform_real_generate_r(0.0, 1.0);
	boost::variate_generator< boost::mt19937&, boost::uniform_real<> > random_real_num_r(generator, uniform_real_generate_r);

	Cluster_Tree_SearchInterval_Set(&(*cur_NRS).Solution_Tree);


	// V5: find out the leaf node to mutate. if the node found is closed then change to another node, 
	// and continue to search for leaf node with dimension =-1

	do{
		min_bound = ((*cur_NRS).Solution_Tree).cur_Interval[(*cur_Node).Dimension][0];
		max_bound = ((*cur_NRS).Solution_Tree).cur_Interval[(*cur_Node).Dimension][1];
		double randr = random_real_num_r();
		//cout << randr << "\t";
		if (randr < 0.5)
		{
			((*cur_NRS).Solution_Tree).cur_Interval[(*cur_Node).Dimension][1] = (*cur_Node).Threshold;
			cur_Node = (*cur_Node).Child_Left;
		}
		else
		{
			((*cur_NRS).Solution_Tree).cur_Interval[(*cur_Node).Dimension][0] = (*cur_Node).Threshold;
			cur_Node = (*cur_Node).Child_Right;
		}
		//cout<< "\n:\t" << ((*cur_NRS).Solution_Tree).cur_Interval[0][0] << "|" << ((*cur_NRS).Solution_Tree).cur_Interval[0][1] << "\t";
		// If the current node is closed, then change to another child node for mutation.
		if ((*cur_Node).Dimension == -2)
		{
			if ((*(*cur_Node).Parent).Child_Left == cur_Node)
			{
				//	cout << "->left->" ;
				((*cur_NRS).Solution_Tree).cur_Interval[(*(*cur_Node).Parent).Dimension][0] = (*(*cur_Node).Parent).Threshold;
				((*cur_NRS).Solution_Tree).cur_Interval[(*(*cur_Node).Parent).Dimension][1] = max_bound;
				cur_Node = (*(*cur_Node).Parent).Child_Right;

			}
			else
			{
				//cout << "->right->";
				((*cur_NRS).Solution_Tree).cur_Interval[(*(*cur_Node).Parent).Dimension][0] = min_bound;
				((*cur_NRS).Solution_Tree).cur_Interval[(*(*cur_Node).Parent).Dimension][1] = (*(*cur_Node).Parent).Threshold;
				cur_Node = (*(*cur_Node).Parent).Child_Left;
			}
		}

		//cout << "\t" << ((*cur_NRS).Solution_Tree).cur_Interval[0][0] << "|" << ((*cur_NRS).Solution_Tree).cur_Interval[0][1] << endl;
	} while ((*cur_Node).Dimension > -1);
	
	for (int j = 0; j < (*cur_NRS).Solution_Tree.no_Dimension; j++)
	{
		(*cur_Chromosome)[j] = (*cur_Node).Optimal_Data[j];
	}
	//printf("before mutate:\n");
	//for (i = 0; i < (*cur_NRS).Solution_Tree.no_Dimension; ++i)
	//{
	//	printf("[%.1f, %.1f]:%.1f\t", (*cur_NRS).Solution_Tree.cur_Interval[i][0], (*cur_NRS).Solution_Tree.cur_Interval[i][1], (*cur_Chromosome)[i]);
	//}
	//printf("\n optimal:\n");
	//for (i = 0; i < (*cur_NRS).Solution_Tree.no_Dimension; ++i)
	//{
	//	printf("%.1f\t", (*cur_Node).Optimal_Data[i]);
	//}

	do {
		//逐维检查是否还剩余可以划分的空间
		for (i = 0; i < (*cur_NRS).Solution_Tree.no_Dimension; ++i)
		{
			if ((*cur_NRS).Solution_Tree.cur_Interval[i][1] - (*cur_NRS).Solution_Tree.cur_Interval[i][0] > 1.0)	 //if there is any space to change																													 //					(*cur_Chromosome)[i] = (int) ((double) rand() / (RAND_MAX + 1) * ((*cur_NRS).Solution_Tree.cur_Interval[i][1] - (*cur_NRS).Solution_Tree.cur_Interval[i][0])) + (*cur_NRS).Solution_Tree.cur_Interval[i][0] + 0.5;																										 //如果有，则变异当前个体。划分当前子空间
				//(*cur_Chromosome)[i] = (int)((double)rand() / (RAND_MAX + 1.0) * ((*cur_NRS).Solution_Tree.cur_Interval[i][1] - (*cur_NRS).Solution_Tree.cur_Interval[i][0])) + (*cur_NRS).Solution_Tree.cur_Interval[i][0] + 0.5;
			{
				double min_bound = (*cur_NRS).Solution_Tree.cur_Interval[i][0] + 0.5;
				double max_bound = (*cur_NRS).Solution_Tree.cur_Interval[i][1] + 0.5;
				boost::uniform_real<> uniform_real_generate_bound(min_bound, max_bound);
				boost::variate_generator< boost::mt19937&, boost::uniform_real<> > random_real_num_bound(generator, uniform_real_generate_bound);
				double rand;
				do {
					rand = random_real_num_bound();
				} while (rand == max_bound);
				(*cur_Chromosome)[i] = (int)rand;
			}
		}
		// 逐维检查是否有维度变异了。如果变异了则跳出循环
		for (i = 0; i < (*cur_NRS).Solution_Tree.no_Dimension; ++i)
		{
			if ((*cur_Chromosome)[i] != (*cur_Node).Optimal_Data[i])  //Optimal_Data = the evaluated point in this partition
				break;												 //once a difference occur, break the mutation
		}
	} while (i == (*cur_NRS).Solution_Tree.no_Dimension);
	//进行检查，变异产生的值是否超出了边界
	for (i = 0; i < (*cur_NRS).Solution_Tree.no_Dimension; ++i)
		if ((*cur_Chromosome)[i] < (*cur_NRS).Solution_Tree.Interval[i][0] || (*cur_Chromosome)[i] > (*cur_NRS).Solution_Tree.Interval[i][1])
		{
			printf("%lf %lf %lf Error !!!\n", (*cur_Chromosome)[i], (*cur_NRS).Solution_Tree.cur_Interval[i][0], (*cur_NRS).Solution_Tree.cur_Interval[i][1]);
			Pause();
		}

	//printf("\n after mutate:\n");
	//for (i = 0; i < (*cur_NRS).Solution_Tree.no_Dimension; ++i)
	//{
	//	printf("[%.1f, %.1f]:%.1f\t", (*cur_NRS).Solution_Tree.cur_Interval[i][0], (*cur_NRS).Solution_Tree.cur_Interval[i][1], (*cur_Chromosome)[i]);
	//}
	//system("pause");
	return;
}
// Input individual
// Output mutate individual
int VisitStrategy_Reintialization(double *cur_Chromosome, double *mutated_Chromosome,
	int visitTimes, int index, NonRevisiting_GA *cur_GA, int axis_range,
	int dim, double *LBound, double *UBound, boost::mt19937 &generator, vector<Selection_Strategy> &Leaf_Vector, const double &mutation_ratio,
	const int Fes, const int maxFes)
{
	int strategy;
	boost::uniform_real<> uniform_real_generate_r(0.0, 1.0);
	boost::variate_generator< boost::mt19937&, boost::uniform_real<> > random_real_num_r(generator, uniform_real_generate_r);
	double tmp_rand = random_real_num_r();
	if ((*(*cur_GA).NRS_Info.Solution_Tree.Root).Dimension == -2)
		//undone: 变异选择的条件
		//mutate to nearest spare space
	{
		Mutate_BSP_FULL++;
		Cal_Archive_Leaf_Probility(Leaf_Vector);
		double tmp_rand = random_real_num_r();
		Cluster_Node* selection_Node = Selection_Roulette(tmp_rand, Leaf_Vector);
		VisitStrategy_On_AxisRange((*selection_Node).Optimal_Data, mutated_Chromosome, dim, generator, axis_range);

		//if (Mutate_On_Axis_Range < 2)
		//{
		//	cout << "rand:\t" << tmp_rand << "\tIndex\t" << (*selection_Node).ID_Self << endl;
		//	for (size_t i = 0; i < Leaf_Vector.size() && i < 20; i++)
		//	{
		//		cout << i << "th: "<< (*Leaf_Vector[i].cur_Node).ID_Self<<"th:\t"<< (*Leaf_Vector[i].cur_Node).Optimal_Fitness<<":\t"
		//			<<Leaf_Vector[i].probability <<"|" << Leaf_Vector[i].accum_probability<<endl;
		//	}
		//	cout << endl;
		//}


		//if (Mutate_On_Axis_Range < 10)
		//{
		//	cout << "chromosome:\t" << (*selection_Node).Optimal_Fitness << "\t";
		//	for (int k = 0; k < dim; k++)
		//		cout << (*selection_Node).Optimal_Data[k] << "\t";
		//	cout << endl;
		//	cout << "mutation:\t";
		//	for (int k = 0; k < dim; k++)
		//		cout << mutated_Chromosome[k] << "\t";
		//	cout << endl;
		//}

		//if (Mutate_On_Axis_Range < 10)
		//{
		//	cout << "Leaf Archive:" << endl;
		//	for (size_t k = 0; k < Leaf_Vector.size(); k++)
		//		cout << (*Leaf_Vector[k].cur_Node).ID_Self << ":\t" << Leaf_Vector[k].probability << "\t" << Leaf_Vector[k].accum_probability << endl;
		//	cout << "selection:\t" << tmp_rand << "\t" << (*selection_Node).ID_Self << endl;
		//}
		strategy = 2;
		Mutate_On_Axis_Range++;
	}
	else
	{
		//		NonRevisitingGA_UnConstraint_Mutation(&(*cur_GA).GA_Info.cur_Population[index], (*cur_GA).GA_Info.no_Dimension, &(*cur_GA).NRS_Info);
		//NonRevisitingGA_UnConstraint_Mutation(&cur_Chromosome, (*cur_GA).GA_Info.no_Dimension, &(*cur_GA).NRS_Info);
		NonRevisitingGA_UnConstraint_Mutation_V4(&cur_Chromosome, (*cur_GA).GA_Info.no_Dimension, &(*cur_GA).NRS_Info, generator);
		for (int j = 0; j < dim; j++)
		{
			mutated_Chromosome[j] = cur_Chromosome[j];
		}
		//		tmp_Chromosome = cur_Chromosome;
		strategy = 1;
		Mutate_To_Nearest_Space++;
	}
	return strategy;
}
void R3PSO_LHC_VisitStrategy(NonRevisiting_GA *cur_GA, int axis_range, double **population, double*population_result, double **speed, double**personal_best,
	double*personal_best_result, int *neighbor_best_index, double *UBound, double *LBound,
	std::vector<std::vector<double> > &archive, std::vector<double> &fitnessArchive, double &maxFitness,
	int population_size, int dim, int &Fes, CEC2013 *pFunc,
	double *mutated_Chromosome, boost::mt19937 &generator, vector<Selection_Strategy> &Leaf_Vector, const double &radius, const double &mutation_ratio)
{
	double dist = 0.0;
	// double radius = pFunc->get_rho();

	int strategy;
	int visitTimes;
	int i, j;
	int best_Archive_Index;
	//double Cur_Chromosome_Fitness;

	for (i = 0; i + 3 <= population_size; i += 3)
	{
		dist = Distance(population[i], population[i + 1], dim);
		dist += Distance(population[i], population[i + 2], dim);
		dist += Distance(population[i + 1], population[i + 2], dim);

		if (sqrt(dist) / 3.0 < radius)
		{
			Converge_Times += 3;
			best_Archive_Index = i;
			best_Archive_Index = personal_best_result[i] > personal_best_result[i + 1] ? i : i + 1;
			best_Archive_Index = personal_best_result[best_Archive_Index] > personal_best_result[i + 2] ? best_Archive_Index : i + 2;
			Archive(personal_best[best_Archive_Index], personal_best_result[best_Archive_Index], archive, fitnessArchive, maxFitness, dim, radius);

			for (j = i; j < i + 3; j++)
			{
				for (int k = 0; k < dim; k++)
				{
					(*cur_GA).GA_Info.cur_Population[j][k] = (int)((double)(population[j][k] - LBound[k]) * (double)(axis_range - 1) / (double)(UBound[k] - LBound[k]));
				}
				visitTimes = VisitStrateg_Record((*cur_GA).GA_Info.cur_Population[j], population_result[j], j,
					cur_GA, &(*cur_GA).NRS_Info.Solution_Tree.Current_ID, (*cur_GA).NRS_Info.Solution_Tree.cur_Interval, Leaf_Vector);

				// Before Mutate
				strategy = VisitStrategy_Reintialization((*cur_GA).GA_Info.cur_Population[j], mutated_Chromosome,
					visitTimes, j, cur_GA, axis_range, dim,
					LBound, UBound, generator, Leaf_Vector, mutation_ratio, Fes, pFunc->get_maxfes());


				for (int k = 0; k < dim; k++)
				{
					population[j][k] = mutated_Chromosome[k] / (double)(axis_range - 1) * (UBound[k] - LBound[k]) + LBound[k];
					//if (population[j][k] > UBound[k])
					//{
					//	cout << population[j][k] << " exceed " << UBound[k] << " chromosom: " << mutated_Chromosome[k] << endl;
					//	system("pause");
					//}
					//if (population[j][k] < LBound[k])
					//{
					//	cout << population[j][k] << " lower " << LBound[k] << " chromosom: " << mutated_Chromosome[k] << endl;
					//	system("pause");
					//}
				}

				memcpy(personal_best[j], population[j], sizeof(double) * dim);
				population_result[j] = Fitness(population[j], Fes, pFunc);

				personal_best_result[j] = population_result[j];

				for (int k = 0; k < dim; k++)
					speed[j][k] = 0.0;
				//After Mutate
				VisitStrateg_Record((*cur_GA).GA_Info.cur_Population[j], population_result[j], j,
					cur_GA, &(*cur_GA).NRS_Info.Solution_Tree.Current_ID, (*cur_GA).NRS_Info.Solution_Tree.cur_Interval, Leaf_Vector);
			}// for(j = i; j < i + 3; j++)
		}// if(sqrt(dist) / 3 < radius)
		else
		{
			for (j = i; j < i + 3; j++)
			{
				//Cur_Chromosome_Fitness = population_result[j];
				for (int k = 0; k < dim; k++)
				{
					(*cur_GA).GA_Info.cur_Population[j][k] = (int)((double)(population[j][k] - LBound[k]) * (double)(axis_range - 1) / (double)(UBound[k] - LBound[k]));
				}

				VisitStrateg_Record((*cur_GA).GA_Info.cur_Population[j], population_result[j], j,
					cur_GA, &(*cur_GA).NRS_Info.Solution_Tree.Current_ID, (*cur_GA).NRS_Info.Solution_Tree.cur_Interval, Leaf_Vector);
			}
		}
	}//	for(i = 0; i+3 <= population_size; i += 3)

}

void R2PSO_LHC_VisitStrategy(NonRevisiting_GA *cur_GA, int axis_range, double **population, double*population_result, double **speed, double**personal_best,
	double*personal_best_result, int *neighbor_best_index, double *UBound, double *LBound,
	std::vector<std::vector<double> > &archive, std::vector<double> &fitnessArchive, double &maxFitness,
	int population_size, int dim, int &Fes, CEC2013 *pFunc,
	double *mutated_Chromosome, boost::mt19937 &generator, vector<Selection_Strategy> &Leaf_Vector, const double &radius, const double &mutation_ratio)
{
	double dist = 0.0;
	int strategy;
	int visitTimes;
	int i, j;
	int best_Archive_Index;
	//double Cur_Chromosome_Fitness;

	for (i = 0; i + 2 <= population_size; i += 2)
	{
		dist = Distance(population[i], population[i + 1], dim);

		if (sqrt(dist) / 2.0 < radius)
		{
			Converge_Times += 2;
			best_Archive_Index = personal_best_result[i] > personal_best_result[i + 1] ? i : i + 1;
			Archive(personal_best[best_Archive_Index], personal_best_result[best_Archive_Index], archive, fitnessArchive, maxFitness, dim, radius);

			for (j = i; j < i + 2; j++)
			{
				for (int k = 0; k < dim; k++)
				{
					(*cur_GA).GA_Info.cur_Population[j][k] = (int)((double)(population[j][k] - LBound[k]) * (double)(axis_range - 1) / (double)(UBound[k] - LBound[k]));
				}
				visitTimes = VisitStrateg_Record((*cur_GA).GA_Info.cur_Population[j], population_result[j], j,
					cur_GA, &(*cur_GA).NRS_Info.Solution_Tree.Current_ID, (*cur_GA).NRS_Info.Solution_Tree.cur_Interval, Leaf_Vector);

				// Before Mutate
				strategy = VisitStrategy_Reintialization((*cur_GA).GA_Info.cur_Population[j], mutated_Chromosome,
					visitTimes, j, cur_GA, axis_range, dim,
					LBound, UBound, generator, Leaf_Vector, mutation_ratio, Fes, pFunc->get_maxfes());
				for (int k = 0; k < dim; k++)
				{
					population[j][k] = mutated_Chromosome[k] / (double)(axis_range - 1) * (UBound[k] - LBound[k]) + LBound[k];
				}
				memcpy(personal_best[j], population[j], sizeof(double) * dim);
				population_result[j] = Fitness(population[j], Fes, pFunc);

				personal_best_result[j] = population_result[j];

				for (int k = 0; k < dim; k++)
					speed[j][k] = 0.0;
				//After Mutate
				VisitStrateg_Record((*cur_GA).GA_Info.cur_Population[j], population_result[j], j,
					cur_GA, &(*cur_GA).NRS_Info.Solution_Tree.Current_ID, (*cur_GA).NRS_Info.Solution_Tree.cur_Interval, Leaf_Vector);
			}// for(j = i; j < i + 2; j++)
		}// if(sqrt(dist) / 2 < radius)
		else
		{
			for (j = i; j < i + 2; j++)
			{
				//Cur_Chromosome_Fitness = population_result[j];
				for (int k = 0; k < dim; k++)
				{
					(*cur_GA).GA_Info.cur_Population[j][k] = (int)((double)(population[j][k] - LBound[k]) * (double)(axis_range - 1) / (double)(UBound[k] - LBound[k]));
				}

				VisitStrateg_Record((*cur_GA).GA_Info.cur_Population[j], population_result[j], j,
					cur_GA, &(*cur_GA).NRS_Info.Solution_Tree.Current_ID, (*cur_GA).NRS_Info.Solution_Tree.cur_Interval, Leaf_Vector);
			}
		}
	}//	for(i = 0; i+2 <= population_size; i += 2)
}

void R3PSO_LHC_VisitStrategy_Random(NonRevisiting_GA *cur_GA, int axis_range, double **population, double*population_result, double **speed, double**personal_best,
	double*personal_best_result, int *neighbor_best_index, double *UBound, double *LBound,
	std::vector<std::vector<double> > &archive, std::vector<double> &fitnessArchive, double &maxFitness,
	int population_size, int dim, int &Fes, CEC2013 *pFunc,
	double *mutated_Chromosome, boost::mt19937 &generator, vector<Selection_Strategy> &Leaf_Vector, const double &radius, const double &mutation_ratio)
{
	double dist = 0.0;


	//int strategy;
	//int visitTimes;
	int i, j;
	int best_Archive_Index;
	//double Cur_Chromosome_Fitness;

	for (i = 0; i + 3 <= population_size; i += 3)
	{
		dist = Distance(population[i], population[i + 1], dim);
		dist += Distance(population[i], population[i + 2], dim);
		dist += Distance(population[i + 1], population[i + 2], dim);

		if (sqrt(dist) / 3.0 < radius)
		{
			Converge_Times += 3;
			//cout << "converge" << i << "\t" << i + 1 << "\t" << i + 2 << endl;
			best_Archive_Index = i;
			best_Archive_Index = personal_best_result[i] > personal_best_result[i + 1] ? i : i + 1;
			best_Archive_Index = personal_best_result[best_Archive_Index] > personal_best_result[i + 2] ? best_Archive_Index : i + 2;
			Archive(personal_best[best_Archive_Index], personal_best_result[best_Archive_Index], archive, fitnessArchive, maxFitness, dim, radius);


			for (j = i; j < i + 3; ++j)
			{
				for (int k = 0; k < dim; ++k)
				{
					boost::uniform_real<> uniform_real_generate_x(LBound[k], UBound[k]);
					boost::variate_generator< boost::mt19937&, boost::uniform_real<> > random_real_num_x(generator, uniform_real_generate_x);
					personal_best[j][k] = population[j][k] = random_real_num_x();
					speed[j][k] = 0.0;
				}
				personal_best_result[j] = population_result[j] = Fitness(population[j], Fes, pFunc);
			}// for(j = i; j < i + 3; j++)
		}// if(sqrt(dist) / 3 < radius)
	}//	for(i = 0; i+3 <= population_size; i += 3)
}


void R2PSO_LHC_VisitStrategy_Random(NonRevisiting_GA *cur_GA, int axis_range, double **population, double*population_result, double **speed, double**personal_best,
	double*personal_best_result, int *neighbor_best_index, double *UBound, double *LBound,
	std::vector<std::vector<double> > &archive, std::vector<double> &fitnessArchive, double &maxFitness,
	int population_size, int dim, int &Fes, CEC2013 *pFunc,
	double *mutated_Chromosome, boost::mt19937 &generator, vector<Selection_Strategy> &Leaf_Vector, const double &radius, const double &mutation_ratio)
{
	double dist = 0.0;
	//int strategy;
	//int visitTimes;
	int i, j;
	int best_Archive_Index;
	//double Cur_Chromosome_Fitness;

	for (i = 0; i + 2 <= population_size; i += 2)
	{
		dist = Distance(population[i], population[i + 1], dim);

		if (sqrt(dist) / 2.0 < radius)
		{
			Converge_Times += 2;
			best_Archive_Index = personal_best_result[i] > personal_best_result[i + 1] ? i : i + 1;
			Archive(personal_best[best_Archive_Index], personal_best_result[best_Archive_Index], archive, fitnessArchive, maxFitness, dim, radius);

			for (j = i; j < i + 2; j++)
			{
				for (int k = 0; k < dim; ++k)
				{
					boost::uniform_real<> uniform_real_generate_x(LBound[k], UBound[k]);
					boost::variate_generator< boost::mt19937&, boost::uniform_real<> > random_real_num_x(generator, uniform_real_generate_x);
					personal_best[j][k] = population[j][k] = random_real_num_x();
					speed[j][k] = 0.0;
				}
				personal_best_result[j] = population_result[j] = Fitness(population[j], Fes, pFunc);
			}// for(j = i; j < i + 2; j++)
		}// if(sqrt(dist) / 2 < radius)
	}//	for(i = 0; i+2 <= population_size; i += 2)
}
double Cal_Mean_Dist(double **population, int begin_index, int dim, int length)
{
	double dist = 0.0;
	for (int i = begin_index; i < begin_index + length - 1; i++)
	{
		for (int j = begin_index + 1; j < begin_index + length; j++)
		{
			dist += Distance(population[i], population[j], dim);
		}
	}
	//return sqrt(dist);
	return dist;

}

