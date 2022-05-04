/**
*  Copyright (c) 2019 SCUT
*  All rights reserved.
*
*  @author      Ting Huang 
*  @date		2019-04-18
*
*  @brief       A probabilistic niching evolutionary computation framework(PNF) is proposed to solve MMOPs.
*
*	@reference
*		T. Huang, Y. -J. Gong, W. -N. Chen, H. Wang and J. Zhang, "A Probabilistic Niching Evolutionary Computation Framework Based on Binary Space Partitioning," in IEEE Transactions on Cybernetics, vol. 52, no. 1, pp. 51-64, Jan. 2022, doi: 10.1109/TCYB.2020.2972907.
*/

#include "./CEC2013/cec2013.h"
#include "header.h"
#include <ctime>
#include <cstdio>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <string>
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/cauchy_distribution.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real.hpp>
#include "BSP/ACC_Math_Kernel/ACC_Math_Kernel.h"
#include "BSP/ClusterTree_Kernel/ClusterTree_Kernel.h"
#include "BSP/GeneticAlgorithm_Kernel/GeneticAlgorithm_Kernel.h"
#include "BSP/NonRevisitingGA_Kernel/NonRevisitingGA_Kernel.h"
#include "BSP/NonRevisitingScheme_Kernel/NonRevisitingScheme_Kernel.h"
#include "BSP/Optimization_TestFunction_Kernel/Optimization_TestFunction_Kernel.h"
#include "BSP/RecordRevisitingScheme_Kernel/RecordRevisitingScheme_Kernel.h"
#include "BSP/VisitStrategy_Kernel/VisitStrategy.h"
#include "BSP/SelectionStrategey_Kernel/Selection_Strategy.h"
int No_Leaf_Initialization;
double DIST_INNER_CLUSTER;
double DIST_INTER_CLUSTER;

using namespace std;
int timesOfRun = 50;
int main(int argc, char *argv[])
{

	int funToRun[] = { 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20 };  //function set
	int population_size_set[] = { 100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100 };

	int i, j, k, fun, run;
	CEC2013 *pFunc = NULL;

	int population_size;
	int funNum = 1;
	int dim;
	int Fes;
	int MAX_FES;
	int temp_global_optima_size;
	int global_optima_num;
	int parameter_NP;
	int record_count = 0;
	int parameter_m = 5;
	int axis_range = 100;
	double parameter_radius = 0;
	double parameter_mutation_ratio = 0;
	double *inner_radius_For_Each_Function = new double[20];
	double *inter_radius_For_Each_Function = new double[20];
	memset(inner_radius_For_Each_Function, 0, sizeof(double) * 20);
	memset(inter_radius_For_Each_Function, 0, sizeof(double) * 20);

	if (argc != 1)
	{
		timesOfRun = atoi(argv[1]);						//axis_range
		parameter_NP = atoi(argv[2]);					//neighbor size
		parameter_m = atoi(argv[3]);					//neighbor size
		axis_range = atoi(argv[4]);						//axis_range
		DIST_INNER_CLUSTER = atof(argv[5]);				//float value
		DIST_INTER_CLUSTER = atof(argv[6]);				//float value, rather than the power
		funNum = argc - 7;								// the first argument is the proc itself
		for (int i = 7; i < argc; i++)
			funToRun[i - 7] = atoi(argv[i]);
	}

	int **runs_global_num = new int*[epsilon_set_size];
	int **runs_used_fitness = new int*[epsilon_set_size];

	bool *epsilon_flag = new bool[epsilon_set_size];

	int record_size = 10;

	int *record_fes = new int[record_size];

	int ***runs_global_num_vs_fes = new int**[epsilon_set_size];

	for (i = 0; i < epsilon_set_size; ++i)
	{
		runs_global_num[i] = new int[timesOfRun];
		runs_used_fitness[i] = new int[timesOfRun];

		runs_global_num_vs_fes[i] = new int*[timesOfRun];

		for (j = 0; j < timesOfRun; ++j)
		{
			runs_global_num_vs_fes[i][j] = new int[record_size];
		}
	}
	for (i = 0; i < 20; i++)
	{
		population_size_set[i] = parameter_NP;
	}

	double weight = 0.7298;
	double phi = 2.05;
	double **PR = new double*[funNum];
	double **SR = new double*[funNum];
	double **avgEval = new double*[funNum];
	for (i = 0; i < funNum; i++)
	{
		PR[i] = new double[epsilon_set_size];
		SR[i] = new double[epsilon_set_size];
		avgEval[i] = new double[epsilon_set_size];
	}
	for (i = 0; i < funNum; i++)
	{
		for (j = 0; j < epsilon_set_size; j++)
		{
			PR[i][j] = 0.0;
			SR[i][j] = 0.0;
			avgEval[i][j] = 0.0;
		}
	}
	char resultFile[200];
	ofstream outFile;
	double *Mutate_Probability_Times = new double[funNum];
	double *Mutate_On_Axis_Range_Times = new double[funNum];
	double *Mutate_To_Nearest_Space_Times = new double[funNum];
	double *Mean_Mutate_BSP_Full = new double[funNum];
	double *Sum_Archive = new double[funNum];
	double *maxLeaf = new double[funNum];
	double *Mean_Max_No_Vist = new double[funNum];
	double *Mean_Converge_Times = new double[funNum];
	double *Better_Times = new double[funNum];
	double *Middle_Times = new double[funNum];
	double *Worse_Times = new double[funNum];

	for (fun = 0; fun < funNum; fun++)
	{
		cout << "Function " << funToRun[fun] << " Begined!" << endl;

		pFunc = new CEC2013(funToRun[fun]);

		dim = pFunc->get_dimension();

		global_optima_num = pFunc->get_no_goptima();

		MAX_FES = pFunc->get_maxfes();
		double *LBound = new double[dim];
		double *UBound = new double[dim];
		for (i = 0; i < dim; ++i)
		{
			LBound[i] = pFunc->get_lbound(i);
			UBound[i] = pFunc->get_ubound(i);
		}

		for (i = 0; i < record_size; ++i)
		{
			record_fes[i] = (int)((double(i + 1)) / record_size * MAX_FES);
		}

		population_size = population_size_set[funToRun[fun] - 1];
		double **population = new double*[population_size];
		double **speed = new double*[population_size];
		double **personal_best = new double*[population_size];

		for (i = 0; i < population_size; ++i)
		{
			population[i] = new double[dim];
			speed[i] = new double[dim];
			personal_best[i] = new double[dim];
		}

		double *results = new double[population_size];
		double *personal_best_results = new double[population_size];
		double *tmp_Chromosome = new double[dim];

		int *neighbor_best_index = new int[population_size];

		Fes = 0;

		boost::mt19937 generator(time(0)*rand());
		boost::uniform_real<> uniform_real_generate_r(0, 1);
		boost::variate_generator< boost::mt19937&, boost::uniform_real<> > random_real_num_r(generator, uniform_real_generate_r);


		NonRevisiting_GA cur_GA;
		double mutation_range = 0.05;
		int max_archive_size = 2 * population_size * NONREVISITINGGA_MAX_GENERATION;
		int activate_mode = 1;
		Optimization_TestFunction_Dimension = dim;
		maxLeaf[fun] = 0;
		Mutate_Probability_Times[fun] = 0;
		Mutate_On_Axis_Range_Times[fun] = 0;
		Mutate_To_Nearest_Space_Times[fun] = 0;
		Mean_Mutate_BSP_Full[fun] = 0;
		Mean_Max_No_Vist[fun] = 0;
		Mean_Converge_Times[fun] = 0;
		Sum_Archive[fun] = 0;
		No_Leaf_Initialization = 0;
		Better_Times[fun] = 0;

		Middle_Times[fun] = 0;
		Worse_Times[fun] = 0;

		// parameter setting
		double MeanBound = 0.0;
		for (int k = 0; k < dim; k++)
		{
			MeanBound += UBound[k] - LBound[k];
		}
		MeanBound /= (double)dim;
		if (DIST_INNER_CLUSTER < 0)
		{
			inner_radius_For_Each_Function[funToRun[fun] - 1] = pow(0.1, -1 * DIST_INNER_CLUSTER);
		}
		else
		{
			inner_radius_For_Each_Function[funToRun[fun] - 1] = pow(0.1, DIST_INNER_CLUSTER) * dim * MeanBound;
		}

		inter_radius_For_Each_Function[funToRun[fun] - 1]  = DIST_INTER_CLUSTER;

		double inner_radius = inner_radius_For_Each_Function[funToRun[fun] - 1];
		double inter_radius = inter_radius_For_Each_Function[funToRun[fun] - 1];
		double mutation_ratio = parameter_mutation_ratio;

		for (run = 0; run < timesOfRun; run++)
		{
			MiddleCounts = WorseCounts = BetterCounts = 0;
			Max_No_Visit = 1;
			Mutate_Probability = 0;
			Mutate_On_Axis_Range = 0;
			Mutate_To_Nearest_Space = 0;
			No_Leaf_Initialization = 0;
			Converge_Times = 0;
			vector<Selection_Strategy> Leaf_Vector;
			// Archive
			vector<vector<double> > archive;
			vector<double> fitnessArchive;
			cout << "Running the " << run << "th times - F"<< funToRun[fun]<< endl;

			for (i = 0; i < epsilon_set_size; ++i)
			{
				epsilon_flag[i] = false;
			}

			Fes = 0;
			record_count = 0;

			// Initialize population
			for (j = 0; j < dim; ++j)
			{
				boost::uniform_real<> uniform_real_generate_x(LBound[j], UBound[j]);
				boost::variate_generator< boost::mt19937&, boost::uniform_real<> > random_real_num_x(generator, uniform_real_generate_x);

				for (i = 0; i < population_size; ++i)
				{
					personal_best[i][j] = population[i][j] = random_real_num_x();
					speed[i][j] = 0.0;
				}
			}
			Popupation_Fitness(population, population_size, Fes, results, pFunc);
			memcpy(personal_best_results, results, sizeof(double) * population_size);

			NonRevisitingGA_Construction(&cur_GA, Optimization_TestFunction_Dimension, population_size, population_size, 0.5, mutation_range, max_archive_size, activate_mode);
			for (i = 0; i < Optimization_TestFunction_Dimension; ++i)
			{
				cur_GA.GA_Info.SearchSpace[i][0] = 0.0 - 0.5;			// for scaled solutions
				cur_GA.GA_Info.SearchSpace[i][1] = axis_range - 0.5;	//
			}
			cur_GA.GA_Info.no_CrossOver_Point = 0;
			cur_GA.NRS_Info.Activate = activate_mode;


			cur_GA.NRS_Info.Max_no_Leaf = 0;
			for (i = 0; i < cur_GA.GA_Info.no_Dimension; ++i)
				for (j = 0; j < 2; ++j)
					cur_GA.NRS_Info.Solution_Tree.Interval[i][j] = cur_GA.GA_Info.SearchSpace[i][j];

			for (i = 0; i < population_size; i++)
			{
				for (j = 0; j < dim; j++)
				{
					cur_GA.GA_Info.cur_Population[i][j] = (int)((double)(population[i][j] - LBound[j]) * (double)(axis_range - 1) / (double)(UBound[j] - LBound[j]));

				}
				VisitStrateg_Record(cur_GA.GA_Info.cur_Population[i], results[i], i,
					&cur_GA, &cur_GA.NRS_Info.Solution_Tree.Current_ID, cur_GA.NRS_Info.Solution_Tree.cur_Interval, Leaf_Vector);
			}
			No_Leaf_Initialization = cur_GA.NRS_Info.Max_no_Leaf;

			double maxFitness;
			while (Fes < MAX_FES)
			{
				Get_Neighbor_Best_m_LHC(personal_best_results, neighbor_best_index, population_size, parameter_m);
				Evolve(population, speed, results, personal_best, personal_best_results, neighbor_best_index, UBound, LBound, population_size, weight, phi, dim, Fes, pFunc);
				VisitStrategy_V8_RPSO(population, results, speed, personal_best, personal_best_results, neighbor_best_index, UBound, LBound,
					population_size, parameter_m, dim, Fes, pFunc, generator, archive, fitnessArchive, maxFitness, tmp_Chromosome,
					Leaf_Vector, &cur_GA, axis_range, inner_radius, inter_radius);
			}
			char fileName_Leaf[200];


			vector<double> indiv(dim);
			for (i = 0; i < population_size; i++)
			{
				for (j = 0; j < dim; j++)
					indiv[j] = personal_best[i][j];
				archive.push_back(indiv);
				fitnessArchive.push_back(personal_best_results[i]);
			}
			double maxFitness_for_archive;
			if (fitnessArchive.empty() == false)
			{
				maxFitness_for_archive = fitnessArchive[0];
				for (i = 0; i < fitnessArchive.size(); i++)
				{
					if (maxFitness_for_archive < fitnessArchive[i])
						maxFitness_for_archive = fitnessArchive[i];
				}
			}
			NewType *fitness_index = new NewType[fitnessArchive.size()];
			for (int i = 0; i < fitnessArchive.size(); i++)
			{
				fitness_index[i].id = i;
				fitness_index[i].data = fitnessArchive[i];
			}
			sort(fitness_index, fitness_index + fitnessArchive.size(), Compare_NewType_Des);

			vector<vector<double> > evaluation_pop_vec;
			vector<double> evaluation_results_vec;
			for (int ep = 0; ep < epsilon_set_size; ++ep)
			{

				Final_Archive_4_Evaluation_V3(archive, fitnessArchive, evaluation_pop_vec, evaluation_results_vec, fitness_index, maxFitness_for_archive, epsilon_set[ep], inter_radius, dim);

				int finalSize = evaluation_results_vec.size();
				double **finalPop = new double*[finalSize];
				for (i = 0; i < finalSize; i++)
				{
					finalPop[i] = new double[dim];
				}
				double *finalResults = new double[finalSize];
				for (i = 0; i < finalSize; i++)
				{
					finalResults[i] = evaluation_results_vec[i];
					for (j = 0; j < dim; j++)
					{
						finalPop[i][j] = evaluation_pop_vec[i][j];
					}
				}

				temp_global_optima_size = Compute_Global_Optima_Found(finalPop, finalResults, finalSize, dim, epsilon_set[ep], pFunc);
				runs_global_num[ep][run] = temp_global_optima_size;
				runs_used_fitness[ep][run] = MAX_FES;
				if (temp_global_optima_size == global_optima_num)
				{
					epsilon_flag[ep] = true;
				}

				for (i = 0; i < finalSize; ++i)
					delete[]finalPop[i];
				delete[]finalPop;
				delete[]finalResults;
				Sum_Archive[fun] += finalSize;
			}
			delete[]fitness_index;
			for (i = 0; i < epsilon_set_size; ++i)
			{
				PR[fun][i] += runs_global_num[i][run];
				if (epsilon_flag[i])
					SR[fun][i] ++;
				avgEval[fun][i] += runs_used_fitness[i][run];
			}
			maxLeaf[fun] += cur_GA.NRS_Info.Max_no_Leaf;

			Mutate_Probability_Times[fun] += Mutate_Probability;
			Mutate_On_Axis_Range_Times[fun] += Mutate_On_Axis_Range;
			Mutate_To_Nearest_Space_Times[fun] += Mutate_To_Nearest_Space;
			Mean_Mutate_BSP_Full[fun] += Mutate_BSP_FULL;
			Mean_Max_No_Vist[fun] += Max_No_Visit;
			Mean_Converge_Times[fun] += Converge_Times;
			
			Better_Times[fun] += BetterCounts;
			Middle_Times[fun] += MiddleCounts;
			Worse_Times[fun] += WorseCounts;

			NonRevisitingGA_Destruction(cur_GA);
		}// end run



		for (i = 0; i < epsilon_set_size; ++i)
		{
			PR[fun][i] /= (double)global_optima_num;
		}

		char fun_name[10];
		char epsilon_counter[10];
		char axis_num[10];
		char NP_num[10];
		char m_num[10];
		char radius_num_inner[10];
		char radius_num_inter[10];
		snprintf(fun_name, 10, "%d", funToRun[fun]);
		snprintf(axis_num, 10, "%d", axis_range);
		snprintf(radius_num_inner, 10, "%f", DIST_INNER_CLUSTER);
		snprintf(radius_num_inter, 10, "%f", DIST_INTER_CLUSTER);
		snprintf(NP_num, 10, "%d", parameter_NP);
		snprintf(m_num, 10, "%d", parameter_m);

		string *Epsilon_Files = new string[epsilon_set_size];
		ofstream *out_Epsilon = new ofstream[epsilon_set_size];

		for (i = 0; i < epsilon_set_size; ++i)
		{
			snprintf(epsilon_counter, 10, "%d", i + 1);
			Epsilon_Files[i] = "./Results/DATA_PNF_PSO/" + string(epsilon_counter) + "/" + "Final_NumOptima_Epsilon_" + string(epsilon_counter)
				+ "_Func_" + string(fun_name) + "_Axis_" + string(axis_num) + "_Rinner_" + string(radius_num_inner) + "_Rinter_" + string(radius_num_inter) +
				"_NP_" + string(NP_num) + "_m_" +
				string(m_num) + ".txt";
		}

		for (i = 0; i < epsilon_set_size; ++i)
		{
			out_Epsilon[i].open(Epsilon_Files[i].c_str());
		}
		for (i = 0; i < epsilon_set_size; ++i)
		{
			for (j = 0; j < timesOfRun; ++j)
			{
				out_Epsilon[i] << runs_global_num[i][j] << "\t" << runs_used_fitness[i][j] << endl;
			}
		}

		for (i = 0; i < epsilon_set_size; ++i)
		{
			out_Epsilon[i].close();
		}

		Better_Times[fun] /= (double)timesOfRun;
		Middle_Times[fun] /= (double)timesOfRun;
		Worse_Times[fun] /= (double)timesOfRun;

		Mutate_Probability_Times[fun] /= (double)timesOfRun;
		Mutate_On_Axis_Range_Times[fun] /= (double)timesOfRun;
		Mutate_To_Nearest_Space_Times[fun] /= (double)timesOfRun;
		Mean_Mutate_BSP_Full[fun] /= (double)timesOfRun;
		Mean_Max_No_Vist[fun] /= (double)timesOfRun;
		Mean_Converge_Times[fun] /= (double)timesOfRun;
		Sum_Archive[fun] /= (double)timesOfRun * (double)epsilon_set_size;
		cout << "Function " << funToRun[fun] << " Finished!" << endl;

		for (i = 0; i < population_size; ++i)
		{
			delete[]population[i];
			delete[]speed[i];
			delete[]personal_best[i];
		}
		delete[]population;
		delete[]speed;
		delete[]personal_best;
		delete[]personal_best_results;
		delete[]tmp_Chromosome;

		delete[]results;
		delete[]LBound;
		delete[]UBound;
		delete[]neighbor_best_index;
		delete[]Epsilon_Files;
		delete[]out_Epsilon;
	}// end func

	sprintf(resultFile, "./Results/DATA_PNF_DE/PNF_PSO_Axis_%d_Rinner_%f_Rinter_%f_NP_%d_m_%d.txt", axis_range, DIST_INNER_CLUSTER, DIST_INTER_CLUSTER, parameter_NP, parameter_m);
	outFile.open(resultFile, ios::out);

	for (j = 0; j < epsilon_set_size; j++)
	{
		cout << epsilon_set[j] << endl;
		outFile << epsilon_set[j] << endl;
		for (i = 0; i < funNum; i++)
		{
			cout << (double)PR[i][j] / timesOfRun << "\t";
			cout << (double)SR[i][j] / timesOfRun << "\t";
			cout << (double)avgEval[i][j] / timesOfRun;
			cout << endl;

			outFile << (double)PR[i][j] / timesOfRun << "\t";
			outFile << (double)SR[i][j] / timesOfRun << "\t";
			outFile << (double)avgEval[i][j] / timesOfRun;
			outFile << endl;
		}
#if defined(__GNUC__) || defined(linux) ||defined(__linux)
		if (funNum < 20)
		{
			for (; i < 20; i++)
			{
				cout << endl;
				outFile << endl;
			}

		}
#endif

		outFile << endl;
		cout << endl;
	}
	outFile << endl;
	cout << endl;
	for (i = 0; i < funNum; i++)
	{
		cout << i + 1 << "th:\t" << endl;
		cout << "Converge Radius:\t" << inner_radius_For_Each_Function[funToRun[i] - 1] << endl;
		cout << "Mean Axis Range:\t" << Mutate_On_Axis_Range_Times[i] << endl;
		cout << "Mean Nearest Space:\t" << Mutate_To_Nearest_Space_Times[i] << endl;
		cout << "Mean Converge Times:\t" << Mean_Converge_Times[i] << endl;
		cout << "Mean Max No Visit:\t" << Mean_Max_No_Vist[i] << endl;
		cout << "Mean Archive:\t" << Sum_Archive[i] << endl;
		cout << "Max Leaf:\t" << maxLeaf[i] << endl;
		cout << "Better Counts;\t" << Better_Times[i] << endl;
		cout << "Middle Counts;\t" << Middle_Times[i] << endl;
		cout << "Worse Counts;\t" << Worse_Times[i] << endl;
		cout << endl;

		outFile << i + 1 << "th:\t" << endl;
		outFile << inner_radius_For_Each_Function[funToRun[i] - 1] << endl;
		outFile << Mutate_On_Axis_Range_Times[i] << endl;
		outFile << Mutate_To_Nearest_Space_Times[i] << endl;
		outFile << Mean_Converge_Times[i] << endl;
		outFile << Mean_Max_No_Vist[i] << endl;
		outFile << Sum_Archive[i] << endl;
		outFile << maxLeaf[i] << endl;
		outFile << Better_Times[i] << endl;
		outFile << Middle_Times[i] << endl;
		outFile << Worse_Times[i] << endl;
		outFile << endl;
	}
	for (int k = 0; k < 20; k++)
		outFile << inner_radius_For_Each_Function[k] << endl;
	outFile << endl << endl;
	for (int k = 0; k < 20; k++)
		cout << inner_radius_For_Each_Function[k] << endl;
	cout << endl << endl;

	for (int k = 0; k < 20; k++)
		outFile << inter_radius_For_Each_Function[k] << endl;
	outFile << endl << endl;
	for (int k = 0; k < 20; k++)
		cout << inter_radius_For_Each_Function[k] << endl;
	cout << endl << endl;

	outFile.close();
	for (i = 0; i < epsilon_set_size; ++i)
	{
		delete[]runs_global_num[i];
		delete[]runs_used_fitness[i];
		for (j = 0; j < timesOfRun; ++j)
		{
			delete[]runs_global_num_vs_fes[i][j];
		}
		delete[]runs_global_num_vs_fes[i];

	}

	delete[]runs_global_num;
	delete[]runs_used_fitness;
	delete[]epsilon_flag;

	delete[]runs_global_num_vs_fes;

	delete[]record_fes;

	for (i = 0; i < funNum; i++)
	{
		delete[]PR[i];
		delete[]SR[i];
		delete[]avgEval[i];
	}
	delete[]PR;
	delete[]SR;
	delete[]avgEval;
	delete[]Mutate_Probability_Times;
	delete[]maxLeaf;
	delete[]Mutate_On_Axis_Range_Times;
	delete[]Mutate_To_Nearest_Space_Times;
	delete[]Mean_Max_No_Vist;
	delete[]Mean_Converge_Times;
	delete[]Sum_Archive;
	delete[]Better_Times;
	delete[]Middle_Times;
	delete[]Worse_Times;
	delete pFunc;
	delete[]inner_radius_For_Each_Function;
	delete[]inter_radius_For_Each_Function;
#ifdef _MSC_VER
	system("pause");
#endif
	return 0;
}