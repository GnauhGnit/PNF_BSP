
#include "BSP/ACC_Math_Kernel/ACC_Math_Kernel.h"
#include "BSP/ClusterTree_Kernel/ClusterTree_Kernel.h"
#include "BSP/GeneticAlgorithm_Kernel/GeneticAlgorithm_Kernel.h"
#include "BSP/NonRevisitingGA_Kernel/NonRevisitingGA_Kernel.h"
#include "BSP/NonRevisitingScheme_Kernel/NonRevisitingScheme_Kernel.h"
#include "BSP/Optimization_TestFunction_Kernel/Optimization_TestFunction_Kernel.h"
#include "header.h"

#include "./BSP/VisitStrategy_Kernel/VisitStrategy.h"


using namespace std;
int BetterCounts;
int MiddleCounts;
int WorseCounts;
// descending
bool Compare_NewType_Des(NewType data1, NewType data2)
{
	return data1.data > data2.data;
}
// ascending
bool Compare_NewType_Asc(NewType data1, NewType data2)
{
	return data1.data < data2.data;
}
bool Compare_NewType(NewType data1, NewType data2)
{
	return data1.data > data2.data;
}


double Fitness(double *particle, int &Fes, CEC2013 *pFunc)
{
	Fes++;
	return pFunc->evaluate(particle);
}


void Popupation_Fitness(double **population, int population_size, int &Fes, double *result, CEC2013 *pFunc)
{
	int i;
	for (i = 0; i < population_size; ++i)
	{
		result[i] = pFunc->evaluate(population[i]);

	}

	Fes += population_size;
}



double Distance(double *vector1, double *vector2, int dim)
{
	double sum = 0;
	for (int i = 0; i < dim; ++i)
	{
		sum += pow(vector1[i] - vector2[i], 2);
	}

	return sqrt(sum);
}
double Distance(vector<double> vector1, vector<double> vector2, int dim)
{
	double sum = 0;
	for (int i = 0; i < dim; ++i)
	{
		sum += pow(vector1[i] - vector2[i], 2);
	}
	return sqrt(sum);
}
double Distance(double *vector1, vector<double> vector2, int dim)
{
	assert(vector2.size() == dim);

	double sum = 0;
	for (int i = 0; i < dim; ++i)
	{
		sum += pow(vector1[i] - vector2[i], 2);
	}

	return sqrt(sum);
}



int Close_Particle(double *child, double **population, int population_size, int dim)
{
	double min_distance = 1e100;

	int i;

	int particle_index = -1;

	double temp_distance;

	for (i = 0; i < population_size; ++i)
	{
		temp_distance = Distance(child, population[i], dim);
		if (min_distance > temp_distance)
		{
			min_distance = temp_distance;
			particle_index = i;
		}
	}


	return particle_index;

}

void Get_Seeds(double **population, double *population_results, int population_size, int dim, vector<double> &seed_fitness, double radius)
{
	int i, j;
	bool found;
	double dist;

	vector< vector<double> > seeds;


	//sort population
	NewType *temp = new NewType[population_size];
	for (i = 0; i < population_size; ++i)
	{
		temp[i].data = population_results[i];
		temp[i].id = i;
	}


	sort(temp, temp + population_size, Compare_NewType_Des);

	seed_fitness.clear();

	for (i = 0; i < population_size; ++i)
	{
		found = false;

		for (j = 0; j < seeds.size(); ++j)
		{
			dist = Distance(population[temp[i].id], seeds[j], dim);
			if (dist <= radius)
			{
				found = true;
				break;
			}
		}

		if (!found)
		{
			vector<double> temp_seed(population[temp[i].id], population[temp[i].id] + dim);
			seeds.push_back(temp_seed);
			seed_fitness.push_back(population_results[temp[i].id]);
		}

	}

	//  for( i = 0;  i < seeds.size(); ++i )
	//  {
	//if (fabs(seed_fitness[i] - 0) < 0.1)
	//{
	//	cout << i << " th:\t" << seed_fitness[i] << "\t";
	//	for (j = 0; j < dim; j++)
	//	{
	//		cout << seeds[i][j] << "\t";
	//	}
	//	cout << endl;
	//}
	//  }

	delete[]temp;
	seeds.clear();
}





int How_Many_Global_Optima(vector<double> seed_fitness, double epsilon, CEC2013 *pFunc)
{

	int  counter = 0;;
	for (unsigned int i = 0; i < seed_fitness.size(); ++i)
	{
		/* |F_seed - F_goptimum| <= accuracy */
		if (fabs(seed_fitness[i] - pFunc->get_fitness_goptima()) <= epsilon)
		{
			++counter;
		}
		/* save time */
		if (counter == pFunc->get_no_goptima())
			break;
	}
	return counter;

}


int Compute_Global_Optima_Found(double **population, double *population_results, int population_size, int dim, double epsilon, CEC2013 *pFunc)
{
	vector<double> seed_fitness;

	Get_Seeds(population, population_results, population_size, dim, seed_fitness, pFunc->get_rho());

	int counter = How_Many_Global_Optima(seed_fitness, epsilon, pFunc);

	seed_fitness.clear();

	return counter;

}


void Evolve(double **population, double**speed, double*population_result, double**personal_best,
	double*personal_best_result, int *neighbor_best_index, double *UBound, double *LBound,
	int population_size, double weight, double phi, int dim, int &Fes, CEC2013 *pFunc)
{
	int i;
	for (i = 0; i < population_size; ++i)
	{
		Generate_Child(population[i], speed[i], personal_best[i], personal_best[neighbor_best_index[i]], UBound, LBound, weight, phi, dim);
		population_result[i] = Fitness(population[i], Fes, pFunc);
		if (population_result[i] > personal_best_result[i])
		{
			personal_best_result[i] = population_result[i];
			memcpy(personal_best[i], population[i], sizeof(double) * dim);
		}
	}
}


void Generate_Child(double *child, double *child_speed, double *exemplar1, double *exemplar2, double *UBound, double *LBound, double weight, double phi, int dim)
{
	boost::mt19937 generator(time(0)*rand());
	boost::uniform_real<> uniform_real_generate_r(0, phi);
	boost::variate_generator< boost::mt19937&, boost::uniform_real<> > random_real_num_r(generator, uniform_real_generate_r);

	int i;

	for (i = 0; i < dim; ++i)
	{
		child_speed[i] = weight*(child_speed[i] + random_real_num_r() * (exemplar1[i] - child[i]) + random_real_num_r() * (exemplar2[i] - child[i]));
		child[i] += child_speed[i];

		if (child[i] < LBound[i])
			child[i] = LBound[i];

		if (child[i] > UBound[i])
			child[i] = UBound[i];
	}

}


void printPop(double **population, double *population_results, int population_size, int dim)
{
	cout << "population detail: " << endl;
	//population_size = population_size > 6 ? 6 : population_size;
	//dim = dim > 3 ? 3 : population_size;
	for (int i = 0; i < population_size; i++)
	{
		cout << i << "th\t" << population_results[i] << ":\t";
		for (int j = 0; j < dim; j++)
		{
			cout << population[i][j] << "\t";
		}
		cout << endl;
	}
}
void printArchive(vector<vector<double> >archive, vector<double> fitnessArchive, int dim)
{
	cout << "archive detail: " << endl;
	if (archive.size() != 0)
		for (int i = 0; i < archive.size(); i++)
		{
			cout << i << "th\t" << fitnessArchive[i] << ":\t";
			for (int j = 0; j < dim; j++)
			{
				cout << archive[i][j] << "\t";
			}
			cout << endl;
		}
}



void Get_Neighbor_Best_m_LHC(double *personal_best_results, int *neighbor_best_index, int population_size, int m)
{
	int i, j;
	int m_best_index;
	int neighbor_size = m;
	for (i = 0; i < population_size; i += neighbor_size)
	{
		m_best_index = i;
		if ((i + m + m) > population_size)
			neighbor_size = population_size - i;
		for (j = i + 1; j < i + neighbor_size; j++)
		{
			m_best_index = personal_best_results[j] > personal_best_results[m_best_index] ? j : m_best_index;
		}
		for (j = i; j < i + neighbor_size; j++)
		{
			neighbor_best_index[j] = m_best_index;
		}
	}
}



//随机生成质心
void Evolve_NrNSDE_V5_random_centroid(double **population, double **population_new, double*results, double *individual_trial,
	double *individual_best, double &individual_best_result, double *UBound, double *LBound,
	int population_size, int neighbor_size, double F, double CR, int dim, int &Fes, CEC2013 *pFunc,
	boost::mt19937 &generator, double **species, double *species_results,
	std::vector<std::vector<double> > &archive, std::vector<double> &fitnessArchive, double &maxFitness,
	double *mutated_Chromosome, vector<Selection_Strategy> &Leaf_Vector,
	NonRevisiting_GA *cur_GA, int axis_range, const double &radius)
{
	int r1, r2, r3;
	int species_num = population_size / neighbor_size;
	vector<NewType> FitnessVec;
	NewType nType;
	vector<NewType> distVec;
	int *SelectionFlag = new int[population_size];
	memset(SelectionFlag, 0, sizeof(int)*population_size);
	Species_Info *speciesInfo = new Species_Info[species_num];
	double *population_new_results = new double[population_size];

	int species_seed_index;
	double species_seed_fitness;
	int min_dist_index;
	int pos;
	int species_current_size = 0;

	vector<int> pop_index(population_size);
	for (int i = 0; i < population_size; i++)
	{
		pop_index[i] = i;
	}
	random_shuffle(pop_index.begin(), pop_index.end());
	// 1. form species

	for (int s = 0; s < species_num; s++)
	{

		speciesInfo[s].beginIndex = s * neighbor_size;

		speciesInfo[s].length = neighbor_size;

		species_current_size += speciesInfo[s].length;
		if (population_size - species_current_size < neighbor_size)
			speciesInfo[s].length += population_size - species_current_size;

		int index;
		for (index = 0; index < population_size; index++)
		{
			if (SelectionFlag[pop_index[index]] == 0)
				break;
		}
		species_seed_index = pop_index[index];
		species_seed_fitness = results[pop_index[index]];
		speciesInfo[s].species_seed_fitness = species_seed_fitness;

		distVec.clear();
		for (int i = 0; i < population_size; i++)
		{
			if (SelectionFlag[pop_index[i]] == 0)
			{
				nType.data = Distance(population[pop_index[i]], population[species_seed_index], dim);
				nType.id = pop_index[i];
				distVec.push_back(nType);
			}
		}
		sort(distVec.begin(), distVec.end(), Compare_NewType_Asc);
		pos = speciesInfo[s].beginIndex;
		speciesInfo[s].index = new int[speciesInfo[s].length];
		for (int i = 0; i < speciesInfo[s].length; i++, pos++)
		{
			memcpy(species[pos], population[distVec[i].id], sizeof(double) * dim);
			species_results[pos] = results[distVec[i].id];
			speciesInfo[s].index[i] = distVec[i].id;
			SelectionFlag[distVec[i].id] = 1;
		}
	}

	for (int i = 0; i < population_size; i++)
	{
		memcpy(population[i], species[i], sizeof(double)*dim);
	}
	memcpy(results, species_results, sizeof(double) * population_size);
	// 2. DE
	boost::uniform_real<> uniform_real_generate_u(0.0, 1.0);
	boost::variate_generator< boost::mt19937&, boost::uniform_real<> > random_real_num_u(generator, uniform_real_generate_u);

	for (int s = 0; s < species_num; s++)
	{
		int species_begin_index = speciesInfo[s].beginIndex;
		int species_length = speciesInfo[s].length;
		for (int i = species_begin_index; i < (species_begin_index + species_length); i++)
		{
			do {
				r1 = (int)(random_real_num_u() * species_length) + species_begin_index;
			} while (r1 == i);
			do {
				r2 = (int)(random_real_num_u() * species_length) + species_begin_index;
			} while (r2 == i || r2 == r1);
			do {
				r3 = (int)(random_real_num_u() * species_length) + species_begin_index;
			} while (r3 == i || r3 == r2 || r3 == r1);

			int n = (int)(random_real_num_u() * dim);
			for (int L = 0; L < dim; L++)
			{
				if (random_real_num_u() < CR || L == (dim - 1))
					individual_trial[n] = population[r1][n] + F * (population[r2][n] - population[r3][n]);
				else
					individual_trial[n] = population[i][n];
				n = (n + 1) % dim;
			}
			// Check boundaries and correct
			for (int j = 0; j < dim; j++)
			{
				if (individual_trial[j] < LBound[j])
					individual_trial[j] = LBound[j];
				else
					if (individual_trial[j] > UBound[j])
						individual_trial[j] = UBound[j];
			}
			double tmpFitness = Fitness(individual_trial, Fes, pFunc);
			memcpy(species[i], individual_trial, sizeof(double)*dim);
			species_results[i] = tmpFitness;
		}	//	for (int i = species_begin_index; i < (species_begin_index + species_length); i++)
	}	//	for (int s = 0; s < species_num; s++)
		// species: new
		// population: old
	for (int s = 0; s < species_num; s++)
	{
		int species_begin_index = speciesInfo[s].beginIndex;
		int species_length = speciesInfo[s].length;
		for (int i = species_begin_index; i < (species_begin_index + species_length); i++)
		{
			if (species_results[i] > results[i])
			{
				memcpy(population[i], species[i], sizeof(double)*dim);
				results[i] = species_results[i];
			}
		}
	}
	delete[]SelectionFlag;
	for (int s = 0; s < species_num; s++)
		delete[]speciesInfo[s].index;
	delete[]speciesInfo;
	delete[]population_new_results;
}
void Evolve_RmPSO_LHC(NonRevisiting_GA *cur_GA, int axis_range, double **population, double*population_result, double **speed, double**personal_best,
	double*personal_best_result, int *neighbor_best_index, double *UBound, double *LBound,
	std::vector<std::vector<double> > &archive, std::vector<double> &fitnessArchive, double &maxFitness,
	int population_size, int dim, int &Fes, CEC2013 *pFunc,
	double *mutated_Chromosome, boost::mt19937 &generator, vector<Selection_Strategy> &Leaf_Vector, const double &radius, const double &mutation_ratio,
	int m)
{
	// m for V8
	double dist = 0.0;
	// double radius = pFunc->get_rho();

	int strategy;
	int visitTimes = 0;
	int i, j;
	int best_Archive_Index;
	//double Cur_Chromosome_Fitness;

	for (i = 0; i + m <= population_size; i += m)
	{
		for (j = i; j < i + m; j++)
		{
			//Cur_Chromosome_Fitness = population_result[j];
			for (int k = 0; k < dim; k++)
			{
				(*cur_GA).GA_Info.cur_Population[j][k] = (int)((double)(population[j][k] - LBound[k]) * (double)(axis_range - 1) / (double)(UBound[k] - LBound[k]));
			}
			VisitStrateg_Record((*cur_GA).GA_Info.cur_Population[j], population_result[j], j,
				cur_GA, &(*cur_GA).NRS_Info.Solution_Tree.Current_ID, (*cur_GA).NRS_Info.Solution_Tree.cur_Interval, Leaf_Vector);
		}

		dist = Cal_Mean_Dist(population, i, dim, m);

		if (dist / (double)m <= radius)
		{
			Converge_Times += m;
			best_Archive_Index = i;
			for (j = i + 1; j < i + m; j++)
			{
				best_Archive_Index = personal_best_result[j] > personal_best_result[best_Archive_Index] ? j : best_Archive_Index;
			}

			Archive(personal_best[best_Archive_Index], personal_best_result[best_Archive_Index], archive, fitnessArchive, maxFitness, dim, radius);

			for (j = i; j < i + m; j++)
			{
				//for (int k = 0; k < dim; k++)
				//{
				//	(*cur_GA).GA_Info.cur_Population[j][k] = (int)((double)(population[j][k] - LBound[k]) * (double)(axis_range - 1) / (double)(UBound[k] - LBound[k]));
				//}
				//visitTimes = VisitStrateg_Record((*cur_GA).GA_Info.cur_Population[j], population_result[j], j,
				//	cur_GA, &(*cur_GA).NRS_Info.Solution_Tree.Current_ID, (*cur_GA).NRS_Info.Solution_Tree.cur_Interval, Leaf_Vector);

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
			}// for(j = i; j < i + 3; j++)
		}// if(sqrt(dist) / 3 < radius)
	}//	for(i = 0; i+3 <= population_size; i += 3)
	if (i != population_size)
	{
		int group_size = population_size - i;
		for (j = i; j < group_size; j++)
		{
			//Cur_Chromosome_Fitness = population_result[j];
			for (int k = 0; k < dim; k++)
			{
				(*cur_GA).GA_Info.cur_Population[j][k] = (int)((double)(population[j][k] - LBound[k]) * (double)(axis_range - 1) / (double)(UBound[k] - LBound[k]));
			}

			VisitStrateg_Record((*cur_GA).GA_Info.cur_Population[j], population_result[j], j,
				cur_GA, &(*cur_GA).NRS_Info.Solution_Tree.Current_ID, (*cur_GA).NRS_Info.Solution_Tree.cur_Interval, Leaf_Vector);
		}
		if (group_size > 1)
		{
			dist = Cal_Mean_Dist(population, i, dim, group_size);
			if (dist / (double)group_size <= radius)
			{
				Converge_Times += group_size;
				best_Archive_Index = i;
				for (j = i + 1; j < i + group_size; j++)
				{
					best_Archive_Index = personal_best_result[j] > personal_best_result[best_Archive_Index] ? j : best_Archive_Index;
				}

				Archive(personal_best[best_Archive_Index], personal_best_result[best_Archive_Index], archive, fitnessArchive, maxFitness, dim, radius);

				for (j = i; j < i + group_size; j++)
				{
					//for (int k = 0; k < dim; k++)
					//{
					//	(*cur_GA).GA_Info.cur_Population[j][k] = (int)((double)(population[j][k] - LBound[k]) * (double)(axis_range - 1) / (double)(UBound[k] - LBound[k]));
					//}
					//visitTimes = VisitStrateg_Record((*cur_GA).GA_Info.cur_Population[j], population_result[j], j,
					//	cur_GA, &(*cur_GA).NRS_Info.Solution_Tree.Current_ID, (*cur_GA).NRS_Info.Solution_Tree.cur_Interval, Leaf_Vector);

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
				}// for(j = i; j < i + 3; j++)
			}// if(sqrt(dist) / 3 < radius)
		}	//if (group_size > 1)
	}
}


void VisitStrategy_V8(double **population, double **population_new, double*population_result, double *individual_trial,
	double *individual_best, double &individual_best_result, double *UBound, double *LBound,
	int population_size, int m, double F, double CR, int dim, int &Fes, CEC2013 *pFunc,
	boost::mt19937 &generator, double **species, double *species_results,
	std::vector<std::vector<double> > &archive, std::vector<double> &fitnessArchive, double &maxFitness,
	double *mutated_Chromosome, vector<Selection_Strategy> &Leaf_Vector,
	NonRevisiting_GA *cur_GA, int axis_range, const double &inner_radius, const double &inter_radius)
{
	double Dist_inner_cluster = inner_radius;
	double Dist_inter_cluster = inter_radius;
	int strategy;
	double radius = 0.0;
	int visitTimes = 0;
	int best_Archive_Index;
	int group_num = population_size / m;
	int randIndex;
	double dist;
	double mean_dist = 0.0;
	//double Cur_Chromosome_Fitness;

	Species_Info2 *species_info = new Species_Info2[group_num];
	double *tmp_individual = new double[dim];
	vector<int> not_convergen_index;

	vector<int> indices_vec;
	for(int i = 0; i < 2*m; i++)
	{
		indices_vec.push_back(i);
	}

	int pos = 0;
	for (int s = 0; s < group_num; s++)
	{
		species_info[s].beginIndex = pos;
		species_info[s].length = m;
		species_info[s].center_pos = new double[dim];
		species_info[s].overlap_index = -1;
		species_info[s].isReInitialize = false;
		memset(species_info[s].center_pos, 0, sizeof(double) * dim);
		if (pos + m + m > population_size)
		{
			species_info[s].length = population_size - pos;
		}

		for (int i = pos; i < pos + species_info[s].length; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				species_info[s].center_pos[j] += population[i][j];
			}

			if (i == pos)
			{
				species_info[s].best_fitness = population_result[pos];
				species_info[s].best_index = pos;
			}
			else
				if (population_result[i] > species_info[s].best_fitness)
				{
					species_info[s].best_fitness = population_result[i];
					species_info[s].best_index = i;
				}
		}	//	for (int i = pos; i < pos + species_info[s].length; i++)
		for (int j = 0; j < dim; j++)
		{
			species_info[s].center_pos[j] /= species_info[s].length;
		}

		pos += species_info[s].length;

		boost::uniform_int<> uniform_int_generate_r(species_info[s].beginIndex, pos - 1);
		boost::variate_generator< boost::mt19937&, boost::uniform_int<> > int_num(generator, uniform_int_generate_r);

		int randIndex = int_num();
		double dist = Distance(species_info[s].center_pos, population[randIndex], dim);
		
		// int conver_select_number = 3;
		// double dist = 0.0;
		// if(conver_select_number > species_info[s].length )
		// 	conver_select_number = species_info[s].length;
		// vector<int> tmp_indices_vec(indices_vec.begin(),indices_vec.begin()+conver_select_number);
		// random_shuffle(tmp_indices_vec.begin(), tmp_indices_vec.end());
		
		// for(int kk =0; kk < tmp_indices_vec.size(); kk++)
		// {
		// 	dist += Distance(species_info[s].center_pos, population[tmp_indices_vec[kk] + species_info[s].beginIndex], dim);
		// }
		// dist/= tmp_indices_vec.size();

		// int conver_select_number = 3;
		// double dist = 0.0;
		// if(conver_select_number > species_info[s].length )
			// conver_select_number = species_info[s].length;
		// if(conver_select_number == species_info[s].length)
		// {
	
			// for (int i = species_info[s].beginIndex; i < species_info[s].length + species_info[s].beginIndex; i++)
			// {
				// dist += Distance(species_info[s].center_pos, population[i], dim);
			// }
			// dist /= m;
		// }else
		// {
			// vector<int> tmp_index;
			// int randIndex;
			// bool exist_flag; 
			// for(int k = 0; k < conver_select_number;)
			// {
				// randIndex = int_num();
				// exist_flag = false;
				// for(int kk =0; kk < tmp_index.size(); kk++)
				// {
					// if(randIndex == tmp_index[kk])
					// {
						// exist_flag = true;
						// break;
					// }
				// }
				// if(exist_flag == false)
				// {
					// tmp_index.push_back(randIndex);
				// }
			// }
				// for(int kk =0; kk < tmp_index.size(); kk++)
				// {
					// dist += Distance(species_info[s].center_pos, population[kk], dim);
				// }
				// dist/= tmp_index.size();

		// }

		//double dist = 0.0;
		//for (int i = species_info[s].beginIndex; i < species_info[s].length + species_info[s].beginIndex; i++)
		//{
		//	dist += Distance(species_info[s].center_pos, population[i], dim);
		//}
		//dist /= m;

		/*cout << "dist:\t" << dist << ":\t" << species_info[s].center_pos[0] << "\t" << population[randIndex][0] << endl;
		cout << population[pos - 5][0] << "\t" << population[pos - 4][0] << "\t" << population[pos - 3][0] << "\t" << population[pos - 2][0] << "\t" << population[pos - 1][0] << endl;*/
		if (dist < Dist_inner_cluster)
		{
			species_info[s].isConvergence = true;
			species_info[s].isReInitialize = true;
		}
		else
		{
			species_info[s].isConvergence = false;
			species_info[s].isReInitialize = false;
			not_convergen_index.push_back(s);
		}

	}	//for (int s = 0; s < group_num; s++)

	if (not_convergen_index.size() <= 2)
	{
#ifdef _MSC_VER
		system("pause");
#endif
		return;
	}
	// calculate min dist 
	int pos_s, pos_ss;
	for (int s = 0; s < not_convergen_index.size(); s++)
	{
		pos_s = not_convergen_index[s];
		species_info[pos_s].min_dist_index = -1;
		for (int ss = 0; ss < not_convergen_index.size(); ss++)
		{
			pos_ss = not_convergen_index[ss];
			if (ss == s)
				continue;
			dist = Distance(species_info[pos_s].center_pos, species_info[pos_ss].center_pos, dim);
			if (species_info[pos_s].min_dist_index == -1 || species_info[pos_s].min_dist > dist)
			{
				species_info[pos_s].min_dist = dist;
				species_info[pos_s].min_dist_index = pos_ss;
			}
		}
	}
	vector<NewType> dist_vec; // group index, min dist
	NewType nt;
	for (int s = 0; s < not_convergen_index.size(); s++)
	{
		nt.id = not_convergen_index[s];
		nt.data = species_info[nt.id].min_dist;
		dist_vec.push_back(nt);
	}
	sort(dist_vec.begin(), dist_vec.end(), Compare_NewType_Asc);

	int better_counts = 0;
	int middle_counts = 0;
	int worse_counts = 0;
	double mutation_ratio = 0.0;
	for (int s = 0; s < dist_vec.size(); s++)
	{
		if (dist_vec[s].data > Dist_inter_cluster)
			break;
		int species_pos = dist_vec[s].id;
		int species_min_pos = species_info[species_pos].min_dist_index;
		species_info[species_pos].overlap_index = species_min_pos;
		if (species_info[species_pos].isReInitialize == true || species_info[species_min_pos].isReInitialize == true)
			continue;
		/*cout << "index:\t" << species_pos << "\t" << species_min_pos << endl;*/


		int max_fitness_index, min_fitness_index;
		if (species_info[species_pos].best_fitness > species_info[species_min_pos].best_fitness)
		{
			max_fitness_index = species_pos;
			min_fitness_index = species_min_pos;
		}
		else
		{
			max_fitness_index = species_min_pos;
			min_fitness_index = species_pos;
		}
		for (int j = 0; j < dim; j++)
		{
			//tmp_individual[j] = (species_info[species_pos].center_pos[j] + species_info[species_min_pos].center_pos[j]) / 2.0;
			tmp_individual[j] = (population[max_fitness_index][j] + population[min_fitness_index][j]) / 2.0;
		}
		double tmp_fitness = Fitness(tmp_individual, Fes, pFunc);


		if (tmp_fitness >= species_info[max_fitness_index].best_fitness)
		{
			pos = species_info[max_fitness_index].beginIndex;
			int worstIndex = pos;
			for (int k = pos + 1; k < pos + species_info[max_fitness_index].length; k++)
			{
				worstIndex = population_result[k] < population_result[worstIndex] ? k : worstIndex;
			}
			memcpy(population[worstIndex], tmp_individual, sizeof(double) * dim);
			population_result[worstIndex] = tmp_fitness;

			species_info[min_fitness_index].isReInitialize = true;
			better_counts++;
			BetterCounts++;
		}
		else if (tmp_fitness < species_info[min_fitness_index].best_fitness)
		{
			worse_counts++;
			WorseCounts++;
			pos = species_info[max_fitness_index].beginIndex;
			int worstIndex = pos;
			for (int k = pos + 1; k < pos + species_info[min_fitness_index].length; k++)
			{
				worstIndex = population_result[k] < population_result[worstIndex] ? k : worstIndex;
			}
			if (tmp_fitness > population_result[worstIndex])
			{
				memcpy(population[worstIndex], tmp_individual, sizeof(double) * dim);
				population_result[worstIndex] = tmp_fitness;
			}
		}
		else
		{
			MiddleCounts++;
			middle_counts++;
			species_info[min_fitness_index].isReInitialize = true;
			pos = species_info[max_fitness_index].beginIndex;
			int worstIndex = pos;
			for (int k = pos + 1; k < pos + species_info[max_fitness_index].length; k++)
			{
				worstIndex = population_result[k] < population_result[worstIndex] ? k : worstIndex;
			}
			if (tmp_fitness > population_result[worstIndex])
			{
				memcpy(population[worstIndex], tmp_individual, sizeof(double) * dim);
				population_result[worstIndex] = tmp_fitness;
			}
		}
	}	// for (int s = 0; s < dist_vec.size(); s++)

		//////record
	for (int i = 0; i < population_size; i++)
	{
		//Cur_Chromosome_Fitness = population_result[j];
		for (int k = 0; k < dim; k++)
		{
			(*cur_GA).GA_Info.cur_Population[i][k] = (int)((double)(population[i][k] - LBound[k]) * (double)(axis_range - 1) / (double)(UBound[k] - LBound[k]));
		}
		VisitStrateg_Record((*cur_GA).GA_Info.cur_Population[i], population_result[i], i,
			cur_GA, &(*cur_GA).NRS_Info.Solution_Tree.Current_ID, (*cur_GA).NRS_Info.Solution_Tree.cur_Interval, Leaf_Vector);
	}

	for (int s = 0; s < group_num; s++)
	{
		if (species_info[s].isReInitialize != true)
			continue;
		pos = species_info[s].beginIndex;

		Converge_Times++;
		best_Archive_Index = pos;
		for (int i = pos + 1; i < pos + species_info[s].length; i++)
		{
			best_Archive_Index = population_result[i] > population_result[best_Archive_Index] ? i : best_Archive_Index;
		}
		Archive(population[best_Archive_Index], population_result[best_Archive_Index], archive, fitnessArchive, maxFitness, dim, radius);


		for (int i = pos; i < pos + species_info[s].length; i++)
		{
			strategy = VisitStrategy_Reintialization((*cur_GA).GA_Info.cur_Population[i], mutated_Chromosome,
				visitTimes, i, cur_GA, axis_range, dim,
				LBound, UBound, generator, Leaf_Vector, mutation_ratio, Fes, pFunc->get_maxfes());
			for (int k = 0; k < dim; k++)
			{
				population[i][k] = mutated_Chromosome[k] / (double)(axis_range - 1) * (UBound[k] - LBound[k]) + LBound[k];
				if (population[i][k] < LBound[k])
					population[i][k] = LBound[k];
				if (population[i][k] > UBound[k])
					population[i][k] = UBound[k];
			}
			population_result[i] = Fitness(population[i], Fes, pFunc);
			//After Mutate
			VisitStrateg_Record((*cur_GA).GA_Info.cur_Population[i], population_result[i], i,
				cur_GA, &(*cur_GA).NRS_Info.Solution_Tree.Current_ID, (*cur_GA).NRS_Info.Solution_Tree.cur_Interval, Leaf_Vector);
		}
	}

	for (int i = 0; i < group_num; i++)
		delete[]species_info[i].center_pos;
	delete[]species_info;
	delete[]tmp_individual;
}

void VisitStrategy_V8_random_reInti(double **population, double **population_new, double*population_result, double *individual_trial,
	double *individual_best, double &individual_best_result, double *UBound, double *LBound,
	int population_size, int m, double F, double CR, int dim, int &Fes, CEC2013 *pFunc,
	boost::mt19937 &generator, double **species, double *species_results,
	std::vector<std::vector<double> > &archive, std::vector<double> &fitnessArchive, double &maxFitness,
	double *mutated_Chromosome, vector<Selection_Strategy> &Leaf_Vector,
	NonRevisiting_GA *cur_GA, int axis_range, const double &inner_radius, const double &inter_radius)
{
	double Dist_inner_cluster = inner_radius;
	double Dist_inter_cluster = inter_radius;
	int strategy;
	double radius = 0.0;
	int visitTimes = 0;
	int best_Archive_Index;
	int group_num = population_size / m;
	int randIndex;
	double dist;
	double mean_dist = 0.0;
	//double Cur_Chromosome_Fitness;

	Species_Info2 *species_info = new Species_Info2[group_num];
	double *tmp_individual = new double[dim];
	vector<int> not_convergen_index;

	int pos = 0;
	for (int s = 0; s < group_num; s++)
	{
		species_info[s].beginIndex = pos;
		species_info[s].length = m;
		species_info[s].center_pos = new double[dim];
		species_info[s].overlap_index = -1;
		species_info[s].isReInitialize = false;
		memset(species_info[s].center_pos, 0, sizeof(double) * dim);
		if (pos + m + m > population_size)
		{
			species_info[s].length = population_size - pos;
		}

		for (int i = pos; i < pos + species_info[s].length; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				species_info[s].center_pos[j] += population[i][j];
			}

			if (i == pos)
			{
				species_info[s].best_fitness = population_result[pos];
				species_info[s].best_index = pos;
			}
			else
				if (population_result[i] > species_info[s].best_fitness)
				{
					species_info[s].best_fitness = population_result[i];
					species_info[s].best_index = i;
				}
		}	//	for (int i = pos; i < pos + species_info[s].length; i++)
		for (int j = 0; j < dim; j++)
		{
			species_info[s].center_pos[j] /= species_info[s].length;
		}

		pos += species_info[s].length;

		boost::uniform_int<> uniform_int_generate_r(species_info[s].beginIndex, pos - 1);
		boost::variate_generator< boost::mt19937&, boost::uniform_int<> > int_num(generator, uniform_int_generate_r);

		int randIndex = int_num();
		double dist = Distance(species_info[s].center_pos, population[randIndex], dim);
		if (dist < Dist_inner_cluster)
		{
			species_info[s].isConvergence = true;
			species_info[s].isReInitialize = true;
		}
		else
		{
			species_info[s].isConvergence = false;
			species_info[s].isReInitialize = false;
			not_convergen_index.push_back(s);
		}

	}	//for (int s = 0; s < group_num; s++)

	if (not_convergen_index.size() <= 2)
	{
#ifdef _MSC_VER
		system("pause");
#endif
		return;
	}
	// calculate min dist 
	int pos_s, pos_ss;
	for (int s = 0; s < not_convergen_index.size(); s++)
	{
		pos_s = not_convergen_index[s];
		species_info[pos_s].min_dist_index = -1;
		for (int ss = 0; ss < not_convergen_index.size(); ss++)
		{
			pos_ss = not_convergen_index[ss];
			if (ss == s)
				continue;
			dist = Distance(species_info[pos_s].center_pos, species_info[pos_ss].center_pos, dim);
			if (species_info[pos_s].min_dist_index == -1 || species_info[pos_s].min_dist > dist)
			{
				species_info[pos_s].min_dist = dist;
				species_info[pos_s].min_dist_index = pos_ss;
			}
		}
	}
	vector<NewType> dist_vec; // group index, min dist
	NewType nt;
	for (int s = 0; s < not_convergen_index.size(); s++)
	{
		nt.id = not_convergen_index[s];
		nt.data = species_info[nt.id].min_dist;
		dist_vec.push_back(nt);
	}
	sort(dist_vec.begin(), dist_vec.end(), Compare_NewType_Asc);

	int better_counts = 0;
	int middle_counts = 0;
	int worse_counts = 0;
	double mutation_ratio = 0.0;
	for (int s = 0; s < dist_vec.size(); s++)
	{
		if (dist_vec[s].data > Dist_inter_cluster)
			break;
		int species_pos = dist_vec[s].id;
		int species_min_pos = species_info[species_pos].min_dist_index;
		species_info[species_pos].overlap_index = species_min_pos;
		if (species_info[species_pos].isReInitialize == true || species_info[species_min_pos].isReInitialize == true)
			continue;
		int max_fitness_index, min_fitness_index;
		if (species_info[species_pos].best_fitness > species_info[species_min_pos].best_fitness)
		{
			max_fitness_index = species_pos;
			min_fitness_index = species_min_pos;
		}
		else
		{
			max_fitness_index = species_min_pos;
			min_fitness_index = species_pos;
		}
		for (int j = 0; j < dim; j++)
		{
			tmp_individual[j] = (population[max_fitness_index][j] + population[min_fitness_index][j]) / 2.0;
		}
		double tmp_fitness = Fitness(tmp_individual, Fes, pFunc);


		if (tmp_fitness >= species_info[max_fitness_index].best_fitness)
		{
			pos = species_info[max_fitness_index].beginIndex;
			int worstIndex = pos;
			for (int k = pos + 1; k < pos + species_info[max_fitness_index].length; k++)
			{
				worstIndex = population_result[k] < population_result[worstIndex] ? k : worstIndex;
			}
			memcpy(population[worstIndex], tmp_individual, sizeof(double) * dim);
			population_result[worstIndex] = tmp_fitness;

			species_info[min_fitness_index].isReInitialize = true;
			better_counts++;
			BetterCounts++;
		}
		else if (tmp_fitness < species_info[min_fitness_index].best_fitness)
		{
			worse_counts++;
			WorseCounts++;
			pos = species_info[max_fitness_index].beginIndex;
			int worstIndex = pos;
			for (int k = pos + 1; k < pos + species_info[min_fitness_index].length; k++)
			{
				worstIndex = population_result[k] < population_result[worstIndex] ? k : worstIndex;
			}
			if (tmp_fitness > population_result[worstIndex])
			{
				memcpy(population[worstIndex], tmp_individual, sizeof(double) * dim);
				population_result[worstIndex] = tmp_fitness;
			}
		}
		else
		{
			MiddleCounts++;
			middle_counts++;
			species_info[min_fitness_index].isReInitialize = true;
			pos = species_info[max_fitness_index].beginIndex;
			int worstIndex = pos;
			for (int k = pos + 1; k < pos + species_info[max_fitness_index].length; k++)
			{
				worstIndex = population_result[k] < population_result[worstIndex] ? k : worstIndex;
			}
			if (tmp_fitness > population_result[worstIndex])
			{
				memcpy(population[worstIndex], tmp_individual, sizeof(double) * dim);
				population_result[worstIndex] = tmp_fitness;
			}
		}
	}	// for (int s = 0; s < dist_vec.size(); s++)

		//////record
	for (int i = 0; i < population_size; i++)
	{
		//Cur_Chromosome_Fitness = population_result[j];
		for (int k = 0; k < dim; k++)
		{
			(*cur_GA).GA_Info.cur_Population[i][k] = (int)((double)(population[i][k] - LBound[k]) * (double)(axis_range - 1) / (double)(UBound[k] - LBound[k]));
		}
		VisitStrateg_Record((*cur_GA).GA_Info.cur_Population[i], population_result[i], i,
			cur_GA, &(*cur_GA).NRS_Info.Solution_Tree.Current_ID, (*cur_GA).NRS_Info.Solution_Tree.cur_Interval, Leaf_Vector);
	}

	for (int s = 0; s < group_num; s++)
	{
		if (species_info[s].isReInitialize != true)
			continue;
		pos = species_info[s].beginIndex;

		Converge_Times++;
		best_Archive_Index = pos;
		for (int i = pos + 1; i < pos + species_info[s].length; i++)
		{
			best_Archive_Index = population_result[i] > population_result[best_Archive_Index] ? i : best_Archive_Index;
		}
		Archive(population[best_Archive_Index], population_result[best_Archive_Index], archive, fitnessArchive, maxFitness, dim, radius);


		for (int i = pos; i < pos + species_info[s].length; i++)
		{

			for (int k = 0; k < dim; k++)
			{
				boost::uniform_real<> uniform_real_generate_x(LBound[k], UBound[k]);
				boost::variate_generator< boost::mt19937&, boost::uniform_real<> > random_real_num_x(generator, uniform_real_generate_x);
				population[i][k] = random_real_num_x();

			}
			population_result[i] = Fitness(population[i], Fes, pFunc);

			//strategy = VisitStrategy_Reintialization((*cur_GA).GA_Info.cur_Population[i], mutated_Chromosome,
			//	visitTimes, i, cur_GA, axis_range, dim,
			//	LBound, UBound, generator, Leaf_Vector, mutation_ratio, Fes, pFunc->get_maxfes());
			//for (int k = 0; k < dim; k++)
			//{
			//	population[i][k] = mutated_Chromosome[k] / (double)(axis_range - 1) * (UBound[k] - LBound[k]) + LBound[k];
			//	if (population[i][k] < LBound[k])
			//		population[i][k] = LBound[k];
			//	if (population[i][k] > UBound[k])
			//		population[i][k] = UBound[k];
			//}
			//population_result[i] = Fitness(population[i], Fes, pFunc);
			//After Mutate
			VisitStrateg_Record((*cur_GA).GA_Info.cur_Population[i], population_result[i], i,
				cur_GA, &(*cur_GA).NRS_Info.Solution_Tree.Current_ID, (*cur_GA).NRS_Info.Solution_Tree.cur_Interval, Leaf_Vector);
		}
	}

	for (int i = 0; i < group_num; i++)
		delete[]species_info[i].center_pos;
	delete[]species_info;
	delete[]tmp_individual;
}

void VisitStrategy_V8_RPSO(double **population, double*population_result,
	double **speed, double**personal_best, double*personal_best_result, int *neighbor_best_index,
	double *UBound, double *LBound, int population_size, int m, int dim, int &Fes, CEC2013 *pFunc,
	boost::mt19937 &generator, std::vector<std::vector<double> > &archive, std::vector<double> &fitnessArchive, double &maxFitness,
	double *mutated_Chromosome, vector<Selection_Strategy> &Leaf_Vector,
	NonRevisiting_GA *cur_GA, int axis_range, const double &inner_radius, const double &inter_radius)
{
	double Dist_inner_cluster = inner_radius;
	double Dist_inter_cluster = inter_radius;
	int strategy;
	double radius = 0.0;
	int visitTimes = 0;
	int best_Archive_Index;
	int group_num = population_size / m;
	int randIndex;
	double dist;
	double mean_dist = 0.0;
	//double Cur_Chromosome_Fitness;

	Species_Info2 *species_info = new Species_Info2[group_num];
	double *tmp_individual = new double[dim];
	vector<int> not_convergen_index;

	int pos = 0;
	for (int s = 0; s < group_num; s++)
	{
		species_info[s].beginIndex = pos;
		species_info[s].length = m;
		species_info[s].center_pos = new double[dim];
		species_info[s].overlap_index = -1;
		species_info[s].isReInitialize = false;
		memset(species_info[s].center_pos, 0, sizeof(double) * dim);
		if (pos + m + m > population_size)
		{
			species_info[s].length = population_size - pos;
		}

		for (int i = pos; i < pos + species_info[s].length; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				species_info[s].center_pos[j] += population[i][j];
			}

			if (i == pos)
			{
				species_info[s].best_fitness = population_result[pos];
				species_info[s].best_index = pos;
			}
			else
				if (population_result[i] > species_info[s].best_fitness)
				{
					species_info[s].best_fitness = population_result[i];
					species_info[s].best_index = i;
				}
		}	//	for (int i = pos; i < pos + species_info[s].length; i++)
		for (int j = 0; j < dim; j++)
		{
			species_info[s].center_pos[j] /= species_info[s].length;
		}

		pos += species_info[s].length;

		boost::uniform_int<> uniform_int_generate_r(species_info[s].beginIndex, pos - 1);
		boost::variate_generator< boost::mt19937&, boost::uniform_int<> > int_num(generator, uniform_int_generate_r);

		int randIndex = int_num();
		double dist = Distance(species_info[s].center_pos, population[randIndex], dim);

		//double dist = 0.0;
		//for (int i = species_info[s].beginIndex; i < species_info[s].length + species_info[s].beginIndex; i++)
		//{
		//	dist += Distance(species_info[s].center_pos, population[i], dim);
		//}
		//dist /= m;

		/*cout << "dist:\t" << dist << ":\t" << species_info[s].center_pos[0] << "\t" << population[randIndex][0] << endl;
		cout << population[pos - 5][0] << "\t" << population[pos - 4][0] << "\t" << population[pos - 3][0] << "\t" << population[pos - 2][0] << "\t" << population[pos - 1][0] << endl;*/
		if (dist < Dist_inner_cluster)
		{
			species_info[s].isConvergence = true;
			species_info[s].isReInitialize = true;
		}
		else
		{
			species_info[s].isConvergence = false;
			species_info[s].isReInitialize = false;
			not_convergen_index.push_back(s);
		}

	}	//for (int s = 0; s < group_num; s++)

	if (not_convergen_index.size() <= 2)
	{
#ifdef _MSC_VER
		system("pause");
#endif
		return;
	}
	// calculate min dist 
	int pos_s, pos_ss;
	for (int s = 0; s < not_convergen_index.size(); s++)
	{
		pos_s = not_convergen_index[s];
		species_info[pos_s].min_dist_index = -1;
		for (int ss = 0; ss < not_convergen_index.size(); ss++)
		{
			pos_ss = not_convergen_index[ss];
			if (ss == s)
				continue;
			dist = Distance(species_info[pos_s].center_pos, species_info[pos_ss].center_pos, dim);
			if (species_info[pos_s].min_dist_index == -1 || species_info[pos_s].min_dist > dist)
			{
				species_info[pos_s].min_dist = dist;
				species_info[pos_s].min_dist_index = pos_ss;
			}
		}
	}
	vector<NewType> dist_vec; // group index, min dist
	NewType nt;
	for (int s = 0; s < not_convergen_index.size(); s++)
	{
		nt.id = not_convergen_index[s];
		nt.data = species_info[nt.id].min_dist;
		dist_vec.push_back(nt);
	}
	sort(dist_vec.begin(), dist_vec.end(), Compare_NewType_Asc);

	int better_counts = 0;
	int middle_counts = 0;
	int worse_counts = 0;
	double mutation_ratio = 0.0;
	for (int s = 0; s < dist_vec.size(); s++)
	{
		if (dist_vec[s].data > Dist_inter_cluster)
			break;
		int species_pos = dist_vec[s].id;
		int species_min_pos = species_info[species_pos].min_dist_index;
		species_info[species_pos].overlap_index = species_min_pos;
		if (species_info[species_pos].isReInitialize == true || species_info[species_min_pos].isReInitialize == true)
			continue;
		/*cout << "index:\t" << species_pos << "\t" << species_min_pos << endl;*/


		int max_fitness_index, min_fitness_index;
		if (species_info[species_pos].best_fitness > species_info[species_min_pos].best_fitness)
		{
			max_fitness_index = species_pos;
			min_fitness_index = species_min_pos;
		}
		else
		{
			max_fitness_index = species_min_pos;
			min_fitness_index = species_pos;
		}
		for (int j = 0; j < dim; j++)
		{
			//tmp_individual[j] = (species_info[species_pos].center_pos[j] + species_info[species_min_pos].center_pos[j]) / 2.0;
			tmp_individual[j] = (population[max_fitness_index][j] + population[min_fitness_index][j]) / 2.0;
		}
		double tmp_fitness = Fitness(tmp_individual, Fes, pFunc);

		if (tmp_fitness >= species_info[max_fitness_index].best_fitness)
		{
			pos = species_info[max_fitness_index].beginIndex;
			int worstIndex = pos;
			for (int k = pos + 1; k < pos + species_info[max_fitness_index].length; k++)
			{
				worstIndex = personal_best_result[k] < personal_best_result[worstIndex] ? k : worstIndex;
			}
			if (tmp_fitness > personal_best_result[worstIndex])
			{
				memcpy(population[worstIndex], tmp_individual, sizeof(double) * dim);
				population_result[worstIndex] = tmp_fitness;

				memcpy(personal_best[worstIndex], tmp_individual, sizeof(double) * dim);
				personal_best_result[worstIndex] = tmp_fitness;
				for (int k = 0; k < dim; k++)
					speed[worstIndex][k] = 0.0;
			}

			species_info[min_fitness_index].isReInitialize = true;
			better_counts++;
			BetterCounts++;
		}
		else if (tmp_fitness < species_info[min_fitness_index].best_fitness)
		{
			worse_counts++;
			WorseCounts++;
			pos = species_info[max_fitness_index].beginIndex;
			int worstIndex = pos;
			for (int k = pos + 1; k < pos + species_info[max_fitness_index].length; k++)
			{
				worstIndex = personal_best_result[k] < personal_best_result[worstIndex] ? k : worstIndex;
			}
			if (tmp_fitness > personal_best_result[worstIndex])
			{
				memcpy(population[worstIndex], tmp_individual, sizeof(double) * dim);
				population_result[worstIndex] = tmp_fitness;

				memcpy(personal_best[worstIndex], tmp_individual, sizeof(double) * dim);
				personal_best_result[worstIndex] = tmp_fitness;
				for (int k = 0; k < dim; k++)
					speed[worstIndex][k] = 0.0;
			}
		}
		else
		{
			MiddleCounts++;
			middle_counts++;
			species_info[min_fitness_index].isReInitialize = true;
			pos = species_info[max_fitness_index].beginIndex;
			int worstIndex = pos;
			for (int k = pos + 1; k < pos + species_info[max_fitness_index].length; k++)
			{
				worstIndex = personal_best_result[k] < personal_best_result[worstIndex] ? k : worstIndex;
			}
			if (tmp_fitness > personal_best_result[worstIndex])
			{
				memcpy(population[worstIndex], tmp_individual, sizeof(double) * dim);
				population_result[worstIndex] = tmp_fitness;

				memcpy(personal_best[worstIndex], tmp_individual, sizeof(double) * dim);
				personal_best_result[worstIndex] = tmp_fitness;
				for (int k = 0; k < dim; k++)
					speed[worstIndex][k] = 0.0;
			}
		}
	}	// for (int s = 0; s < dist_vec.size(); s++)

		//////record
	for (int i = 0; i < population_size; i++)
	{
		//Cur_Chromosome_Fitness = population_result[j];
		for (int k = 0; k < dim; k++)
		{
			(*cur_GA).GA_Info.cur_Population[i][k] = (int)((double)(population[i][k] - LBound[k]) * (double)(axis_range - 1) / (double)(UBound[k] - LBound[k]));
		}
		VisitStrateg_Record((*cur_GA).GA_Info.cur_Population[i], population_result[i], i,
			cur_GA, &(*cur_GA).NRS_Info.Solution_Tree.Current_ID, (*cur_GA).NRS_Info.Solution_Tree.cur_Interval, Leaf_Vector);
	}

	for (int s = 0; s < group_num; s++)
	{
		if (species_info[s].isReInitialize != true)
			continue;
		pos = species_info[s].beginIndex;

		Converge_Times++;
		best_Archive_Index = pos;
		for (int i = pos + 1; i < pos + species_info[s].length; i++)
		{
			best_Archive_Index = personal_best_result[i] > personal_best_result[best_Archive_Index] ? i : best_Archive_Index;
		}
		Archive(personal_best[best_Archive_Index], personal_best_result[best_Archive_Index], archive, fitnessArchive, maxFitness, dim, radius);


		for (int i = pos; i < pos + species_info[s].length; i++)
		{
			strategy = VisitStrategy_Reintialization((*cur_GA).GA_Info.cur_Population[i], mutated_Chromosome,
				visitTimes, i, cur_GA, axis_range, dim,
				LBound, UBound, generator, Leaf_Vector, mutation_ratio, Fes, pFunc->get_maxfes());
			for (int k = 0; k < dim; k++)
			{
				population[i][k] = mutated_Chromosome[k] / (double)(axis_range - 1) * (UBound[k] - LBound[k]) + LBound[k];
				if (population[i][k] < LBound[k])
					population[i][k] = LBound[k];
				if (population[i][k] > UBound[k])
					population[i][k] = UBound[k];
			}
			memcpy(personal_best[i], population[i], sizeof(double) * dim);
			population_result[i] = Fitness(population[i], Fes, pFunc);
			personal_best_result[i] = population_result[i];
			for (int k = 0; k < dim; k++)
				speed[i][k] = 0.0;
			//After Mutate
			VisitStrateg_Record((*cur_GA).GA_Info.cur_Population[i], population_result[i], i,
				cur_GA, &(*cur_GA).NRS_Info.Solution_Tree.Current_ID, (*cur_GA).NRS_Info.Solution_Tree.cur_Interval, Leaf_Vector);
		}
	}
	for (int i = 0; i < group_num; i++)
		delete[]species_info[i].center_pos;
	delete[]species_info;
	delete[]tmp_individual;
}

void VisitStrategy_V8_RPSO_random_reInti(double **population, double*population_result,
	double **speed, double**personal_best, double*personal_best_result, int *neighbor_best_index,
	double *UBound, double *LBound, int population_size, int m, int dim, int &Fes, CEC2013 *pFunc,
	boost::mt19937 &generator, std::vector<std::vector<double> > &archive, std::vector<double> &fitnessArchive, double &maxFitness,
	double *mutated_Chromosome, vector<Selection_Strategy> &Leaf_Vector,
	NonRevisiting_GA *cur_GA, int axis_range, const double &inner_radius, const double &inter_radius)
{
	double Dist_inner_cluster = inner_radius;
	double Dist_inter_cluster = inter_radius;
	int strategy;
	double radius = 0.0;
	int visitTimes = 0;
	int best_Archive_Index;
	int group_num = population_size / m;
	int randIndex;
	double dist;
	double mean_dist = 0.0;
	//double Cur_Chromosome_Fitness;

	Species_Info2 *species_info = new Species_Info2[group_num];
	double *tmp_individual = new double[dim];
	vector<int> not_convergen_index;

	int pos = 0;
	for (int s = 0; s < group_num; s++)
	{
		species_info[s].beginIndex = pos;
		species_info[s].length = m;
		species_info[s].center_pos = new double[dim];
		species_info[s].overlap_index = -1;
		species_info[s].isReInitialize = false;
		memset(species_info[s].center_pos, 0, sizeof(double) * dim);
		if (pos + m + m > population_size)
		{
			species_info[s].length = population_size - pos;
		}

		for (int i = pos; i < pos + species_info[s].length; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				species_info[s].center_pos[j] += population[i][j];
			}

			if (i == pos)
			{
				species_info[s].best_fitness = population_result[pos];
				species_info[s].best_index = pos;
			}
			else
				if (population_result[i] > species_info[s].best_fitness)
				{
					species_info[s].best_fitness = population_result[i];
					species_info[s].best_index = i;
				}
		}	//	for (int i = pos; i < pos + species_info[s].length; i++)
		for (int j = 0; j < dim; j++)
		{
			species_info[s].center_pos[j] /= species_info[s].length;
		}

		pos += species_info[s].length;

		boost::uniform_int<> uniform_int_generate_r(species_info[s].beginIndex, pos - 1);
		boost::variate_generator< boost::mt19937&, boost::uniform_int<> > int_num(generator, uniform_int_generate_r);

		int randIndex = int_num();
		double dist = Distance(species_info[s].center_pos, population[randIndex], dim);

		//double dist = 0.0;
		//for (int i = species_info[s].beginIndex; i < species_info[s].length + species_info[s].beginIndex; i++)
		//{
		//	dist += Distance(species_info[s].center_pos, population[i], dim);
		//}
		//dist /= m;

		/*cout << "dist:\t" << dist << ":\t" << species_info[s].center_pos[0] << "\t" << population[randIndex][0] << endl;
		cout << population[pos - 5][0] << "\t" << population[pos - 4][0] << "\t" << population[pos - 3][0] << "\t" << population[pos - 2][0] << "\t" << population[pos - 1][0] << endl;*/
		if (dist < Dist_inner_cluster)
		{
			species_info[s].isConvergence = true;
			species_info[s].isReInitialize = true;
		}
		else
		{
			species_info[s].isConvergence = false;
			species_info[s].isReInitialize = false;
			not_convergen_index.push_back(s);
		}

	}	//for (int s = 0; s < group_num; s++)

	if (not_convergen_index.size() <= 2)
	{
#ifdef _MSC_VER
		system("pause");
#endif
		return;
	}
	// calculate min dist 
	int pos_s, pos_ss;
	for (int s = 0; s < not_convergen_index.size(); s++)
	{
		pos_s = not_convergen_index[s];
		species_info[pos_s].min_dist_index = -1;
		for (int ss = 0; ss < not_convergen_index.size(); ss++)
		{
			pos_ss = not_convergen_index[ss];
			if (ss == s)
				continue;
			dist = Distance(species_info[pos_s].center_pos, species_info[pos_ss].center_pos, dim);
			if (species_info[pos_s].min_dist_index == -1 || species_info[pos_s].min_dist > dist)
			{
				species_info[pos_s].min_dist = dist;
				species_info[pos_s].min_dist_index = pos_ss;
			}
		}
	}
	vector<NewType> dist_vec; // group index, min dist
	NewType nt;
	for (int s = 0; s < not_convergen_index.size(); s++)
	{
		nt.id = not_convergen_index[s];
		nt.data = species_info[nt.id].min_dist;
		dist_vec.push_back(nt);
	}
	sort(dist_vec.begin(), dist_vec.end(), Compare_NewType_Asc);

	int better_counts = 0;
	int middle_counts = 0;
	int worse_counts = 0;
	double mutation_ratio = 0.0;
	for (int s = 0; s < dist_vec.size(); s++)
	{
		if (dist_vec[s].data > Dist_inter_cluster)
			break;
		int species_pos = dist_vec[s].id;
		int species_min_pos = species_info[species_pos].min_dist_index;
		species_info[species_pos].overlap_index = species_min_pos;
		if (species_info[species_pos].isReInitialize == true || species_info[species_min_pos].isReInitialize == true)
			continue;
		/*cout << "index:\t" << species_pos << "\t" << species_min_pos << endl;*/


		int max_fitness_index, min_fitness_index;
		if (species_info[species_pos].best_fitness > species_info[species_min_pos].best_fitness)
		{
			max_fitness_index = species_pos;
			min_fitness_index = species_min_pos;
		}
		else
		{
			max_fitness_index = species_min_pos;
			min_fitness_index = species_pos;
		}
		for (int j = 0; j < dim; j++)
		{
			//tmp_individual[j] = (species_info[species_pos].center_pos[j] + species_info[species_min_pos].center_pos[j]) / 2.0;
			tmp_individual[j] = (population[max_fitness_index][j] + population[min_fitness_index][j]) / 2.0;
		}
		double tmp_fitness = Fitness(tmp_individual, Fes, pFunc);

		if (tmp_fitness >= species_info[max_fitness_index].best_fitness)
		{
			pos = species_info[max_fitness_index].beginIndex;
			int worstIndex = pos;
			for (int k = pos + 1; k < pos + species_info[max_fitness_index].length; k++)
			{
				worstIndex = personal_best_result[k] < personal_best_result[worstIndex] ? k : worstIndex;
			}
			if (tmp_fitness > personal_best_result[worstIndex])
			{
				memcpy(population[worstIndex], tmp_individual, sizeof(double) * dim);
				population_result[worstIndex] = tmp_fitness;

				memcpy(personal_best[worstIndex], tmp_individual, sizeof(double) * dim);
				personal_best_result[worstIndex] = tmp_fitness;
				for (int k = 0; k < dim; k++)
					speed[worstIndex][k] = 0.0;
			}

			species_info[min_fitness_index].isReInitialize = true;
			better_counts++;
			BetterCounts++;
		}
		else if (tmp_fitness < species_info[min_fitness_index].best_fitness)
		{
			worse_counts++;
			WorseCounts++;
			pos = species_info[max_fitness_index].beginIndex;
			int worstIndex = pos;
			for (int k = pos + 1; k < pos + species_info[max_fitness_index].length; k++)
			{
				worstIndex = personal_best_result[k] < personal_best_result[worstIndex] ? k : worstIndex;
			}
			if (tmp_fitness > personal_best_result[worstIndex])
			{
				memcpy(population[worstIndex], tmp_individual, sizeof(double) * dim);
				population_result[worstIndex] = tmp_fitness;

				memcpy(personal_best[worstIndex], tmp_individual, sizeof(double) * dim);
				personal_best_result[worstIndex] = tmp_fitness;
				for (int k = 0; k < dim; k++)
					speed[worstIndex][k] = 0.0;
			}
		}
		else
		{
			MiddleCounts++;
			middle_counts++;
			species_info[min_fitness_index].isReInitialize = true;
			pos = species_info[max_fitness_index].beginIndex;
			int worstIndex = pos;
			for (int k = pos + 1; k < pos + species_info[max_fitness_index].length; k++)
			{
				worstIndex = personal_best_result[k] < personal_best_result[worstIndex] ? k : worstIndex;
			}
			if (tmp_fitness > personal_best_result[worstIndex])
			{
				memcpy(population[worstIndex], tmp_individual, sizeof(double) * dim);
				population_result[worstIndex] = tmp_fitness;

				memcpy(personal_best[worstIndex], tmp_individual, sizeof(double) * dim);
				personal_best_result[worstIndex] = tmp_fitness;
				for (int k = 0; k < dim; k++)
					speed[worstIndex][k] = 0.0;
			}
		}
	}	// for (int s = 0; s < dist_vec.size(); s++)

		//////record
	for (int i = 0; i < population_size; i++)
	{
		//Cur_Chromosome_Fitness = population_result[j];
		for (int k = 0; k < dim; k++)
		{
			(*cur_GA).GA_Info.cur_Population[i][k] = (int)((double)(population[i][k] - LBound[k]) * (double)(axis_range - 1) / (double)(UBound[k] - LBound[k]));
		}
		VisitStrateg_Record((*cur_GA).GA_Info.cur_Population[i], population_result[i], i,
			cur_GA, &(*cur_GA).NRS_Info.Solution_Tree.Current_ID, (*cur_GA).NRS_Info.Solution_Tree.cur_Interval, Leaf_Vector);
	}

	for (int s = 0; s < group_num; s++)
	{
		if (species_info[s].isReInitialize != true)
			continue;
		pos = species_info[s].beginIndex;

		Converge_Times++;
		best_Archive_Index = pos;
		for (int i = pos + 1; i < pos + species_info[s].length; i++)
		{
			best_Archive_Index = personal_best_result[i] > personal_best_result[best_Archive_Index] ? i : best_Archive_Index;
		}
		Archive(personal_best[best_Archive_Index], personal_best_result[best_Archive_Index], archive, fitnessArchive, maxFitness, dim, radius);


		for (int i = pos; i < pos + species_info[s].length; i++)
		{

			for (int k = 0; k < dim; k++)
			{
				boost::uniform_real<> uniform_real_generate_x(LBound[k], UBound[k]);
				boost::variate_generator< boost::mt19937&, boost::uniform_real<> > random_real_num_x(generator, uniform_real_generate_x);
				population[i][k] = random_real_num_x();
			}
			population_result[i] = Fitness(population[i], Fes, pFunc);

			//strategy = VisitStrategy_Reintialization((*cur_GA).GA_Info.cur_Population[i], mutated_Chromosome,
			//	visitTimes, i, cur_GA, axis_range, dim,
			//	LBound, UBound, generator, Leaf_Vector, mutation_ratio, Fes, pFunc->get_maxfes());
			//for (int k = 0; k < dim; k++)
			//{
			//	population[i][k] = mutated_Chromosome[k] / (double)(axis_range - 1) * (UBound[k] - LBound[k]) + LBound[k];
			//	if (population[i][k] < LBound[k])
			//		population[i][k] = LBound[k];
			//	if (population[i][k] > UBound[k])
			//		population[i][k] = UBound[k];
			//}
			memcpy(personal_best[i], population[i], sizeof(double) * dim);
			personal_best_result[i] = population_result[i];
			for (int k = 0; k < dim; k++)
				speed[i][k] = 0.0;
			//After Mutate
			VisitStrateg_Record((*cur_GA).GA_Info.cur_Population[i], population_result[i], i,
				cur_GA, &(*cur_GA).NRS_Info.Solution_Tree.Current_ID, (*cur_GA).NRS_Info.Solution_Tree.cur_Interval, Leaf_Vector);
		}
	}
	for (int i = 0; i < group_num; i++)
		delete[]species_info[i].center_pos;
	delete[]species_info;
	delete[]tmp_individual;
}


void VisitStrategy_V9(double **population, double **population_new, double*population_result, double *individual_trial,
	double *individual_best, double &individual_best_result, double *UBound, double *LBound,
	int population_size, int m, double F, double CR, int dim, int &Fes, CEC2013 *pFunc,
	boost::mt19937 &generator, double **species, double *species_results,
	std::vector<std::vector<double> > &archive, std::vector<double> &fitnessArchive, double &maxFitness,
	double *mutated_Chromosome, vector<Selection_Strategy> &Leaf_Vector,
	NonRevisiting_GA *cur_GA, int axis_range, const double &inner_radius, const double &inter_radius)
{
	double Dist_inner_cluster = inner_radius;
	double Dist_inter_cluster = inter_radius;
	int strategy;
	double radius = 0.0;
	int visitTimes = 0;
	int best_Archive_Index;
	int group_num = population_size / m;
	int randIndex;
	double dist;
	double mean_dist = 0.0;
	//double Cur_Chromosome_Fitness;

	Species_Info2 *species_info = new Species_Info2[group_num];
	double *tmp_individual = new double[dim];
	vector<int> not_convergen_index;

	int pos = 0;
	for (int s = 0; s < group_num; s++)
	{
		species_info[s].beginIndex = pos;
		species_info[s].length = m;
		species_info[s].center_pos = new double[dim];
		species_info[s].overlap_index = -1;
		species_info[s].isReInitialize = false;
		memset(species_info[s].center_pos, 0, sizeof(double) * dim);
		if (pos + m + m > population_size)
		{
			species_info[s].length = population_size - pos;
		}

		for (int i = pos; i < pos + species_info[s].length; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				species_info[s].center_pos[j] += population[i][j];
			}

			if (i == pos || population_result[i] > species_info[s].best_fitness)
			{
				species_info[s].best_fitness = population_result[i];
				species_info[s].best_index = i;
			}
		}	//	for (int i = pos; i < pos + species_info[s].length; i++)
		for (int j = 0; j < dim; j++)
		{
			species_info[s].center_pos[j] /= species_info[s].length;
		}

		pos += species_info[s].length;

		boost::uniform_int<> uniform_int_generate_r(species_info[s].beginIndex, pos - 1);
		boost::variate_generator< boost::mt19937&, boost::uniform_int<> > int_num(generator, uniform_int_generate_r);

		int randIndex = int_num();
		double dist = Distance(species_info[s].center_pos, population[randIndex], dim);

		if (dist < Dist_inner_cluster)
		{
			species_info[s].isConvergence = true;
			species_info[s].isReInitialize = true;
		}
		else
		{
			species_info[s].isConvergence = false;
			species_info[s].isReInitialize = false;
			//			not_convergen_index.push_back(s);
		}
	}	//for (int s = 0; s < group_num; s++)

	vector<NewType> fitness_vec; // group index, min dist
	NewType nt;
	for (int s = 0; s < group_num; s++)
	{
		nt.id = s;
		nt.data = species_info[s].best_fitness;
		fitness_vec.push_back(nt);
	}
	sort(fitness_vec.begin(), fitness_vec.end(), Compare_NewType_Des);


	// calculate min dist
	int better_counts, middle_counts, worse_counts;
	better_counts = middle_counts = worse_counts = 0;
	int pos_seed, pos_ind;
	for (int s = 0; s < group_num - 1; s++)
	{
		pos_seed = fitness_vec[s].id;
		for (int ss = s + 1; ss < group_num; ss++)
		{
			pos_ind = fitness_vec[ss].id;
			if (species_info[pos_ind].isReInitialize == true)
				continue;
			
			dist = Distance(species_info[pos_seed].center_pos, species_info[pos_ind].center_pos, dim);
			if (dist > Dist_inter_cluster)
				continue;
			for (int j = 0; j < dim; j++)
			{
				tmp_individual[j] = (population[species_info[pos_seed].best_index][j] + population[species_info[pos_ind].best_index][j]) / 2.0;
			}
			double tmp_fitness = Fitness(tmp_individual, Fes, pFunc);

			if (tmp_fitness >= species_info[pos_seed].best_fitness)
			{

				int worstIndex = pos = species_info[pos_seed].beginIndex;
				for (int k = pos + 1; k < pos + species_info[pos_seed].length; k++)
				{
					worstIndex = population_result[k] < population_result[worstIndex] ? k : worstIndex;
				}
				memcpy(population[worstIndex], tmp_individual, sizeof(double) * dim);
				population_result[worstIndex] = tmp_fitness;

				species_info[pos_seed].best_fitness = tmp_fitness;
				species_info[pos_seed].best_index = worstIndex;

				species_info[pos_ind].isReInitialize = true;

				better_counts++;
				BetterCounts++;
			}
			else if (tmp_fitness < species_info[pos_ind].best_fitness)
			{
				int worstIndex = pos = species_info[pos_ind].beginIndex;
				for (int k = pos + 1; k < pos + species_info[pos_ind].length; k++)
				{
					worstIndex = population_result[k] < population_result[worstIndex] ? k : worstIndex;
				}
				if (tmp_fitness > population_result[worstIndex])
				{
					memcpy(population[worstIndex], tmp_individual, sizeof(double) * dim);
					population_result[worstIndex] = tmp_fitness;
				}
				worse_counts++;
				WorseCounts++;
			}
			else
			{
				int worstIndex = pos = species_info[pos_seed].beginIndex;
				for (int k = pos + 1; k < pos + species_info[pos_seed].length; k++)
				{
					worstIndex = population_result[k] < population_result[worstIndex] ? k : worstIndex;
				}
				memcpy(population[worstIndex], tmp_individual, sizeof(double) * dim);
				population_result[worstIndex] = tmp_fitness;

				species_info[pos_seed].best_fitness = tmp_fitness;
				species_info[pos_seed].best_index = worstIndex;

				species_info[pos_ind].isReInitialize = true;

				MiddleCounts++;
				middle_counts++;
			}
		}
	}
	double mutation_ratio;
	//////record
	for (int i = 0; i < population_size; i++)
	{
		//Cur_Chromosome_Fitness = population_result[j];
		for (int k = 0; k < dim; k++)
		{
			(*cur_GA).GA_Info.cur_Population[i][k] = (int)((double)(population[i][k] - LBound[k]) * (double)(axis_range - 1) / (double)(UBound[k] - LBound[k]));
		}
		VisitStrateg_Record((*cur_GA).GA_Info.cur_Population[i], population_result[i], i,
			cur_GA, &(*cur_GA).NRS_Info.Solution_Tree.Current_ID, (*cur_GA).NRS_Info.Solution_Tree.cur_Interval, Leaf_Vector);
	}

	for (int s = 0; s < group_num; s++)
	{
		if (species_info[s].isReInitialize != true)
			continue;
		pos = species_info[s].beginIndex;

		Converge_Times++;
		best_Archive_Index = pos;
		for (int i = pos + 1; i < pos + species_info[s].length; i++)
		{
			best_Archive_Index = population_result[i] > population_result[best_Archive_Index] ? i : best_Archive_Index;
		}
		Archive(population[best_Archive_Index], population_result[best_Archive_Index], archive, fitnessArchive, maxFitness, dim, radius);


		for (int i = pos; i < pos + species_info[s].length; i++)
		{
			strategy = VisitStrategy_Reintialization((*cur_GA).GA_Info.cur_Population[i], mutated_Chromosome,
				visitTimes, i, cur_GA, axis_range, dim,
				LBound, UBound, generator, Leaf_Vector, mutation_ratio, Fes, pFunc->get_maxfes());
			for (int k = 0; k < dim; k++)
			{
				population[i][k] = mutated_Chromosome[k] / (double)(axis_range - 1) * (UBound[k] - LBound[k]) + LBound[k];
				if (population[i][k] < LBound[k])
					population[i][k] = LBound[k];
				if (population[i][k] > UBound[k])
					population[i][k] = UBound[k];
			}
			population_result[i] = Fitness(population[i], Fes, pFunc);
			//After Mutate
			VisitStrateg_Record((*cur_GA).GA_Info.cur_Population[i], population_result[i], i,
				cur_GA, &(*cur_GA).NRS_Info.Solution_Tree.Current_ID, (*cur_GA).NRS_Info.Solution_Tree.cur_Interval, Leaf_Vector);
		}
	}

	for (int i = 0; i < group_num; i++)
		delete[]species_info[i].center_pos;
	delete[]species_info;
	delete[]tmp_individual;
}
