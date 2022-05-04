#ifndef SELF_DEFINE_FUNCTIONS_H_INCLUDED
#define SELF_DEFINE_FUNCTIONS_H_INCLUDED


#include "./CEC2013/cec2013.h"
#include <math.h>
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/cauchy_distribution.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real.hpp>
#include <vector>


#include "./Archive/Archive.h"
#include "./BSP/ACC_Math_Kernel/ACC_Math_Kernel.h" 
#include "./BSP/GeneticAlgorithm_Kernel/GeneticAlgorithm_Kernel.h"
#include "./BSP/NonRevisitingGA_Kernel/NonRevisitingGA_Kernel.h"
#include "./BSP/Optimization_TestFunction_Kernel/Optimization_TestFunction_Kernel.h"

#include "./BSP/RecordRevisitingScheme_Kernel/RecordRevisitingScheme_Kernel.h"
#include "./BSP/SelectionStrategey_Kernel/Selection_Strategy.h"

//#include <ctime>
//#include <cstddef>
using namespace std;

struct Species_Info2
{
	int beginIndex;
	int length;

	double best_fitness;
	int best_index;

	double min_dist;
	int min_dist_index;

	double *center_pos;

	int overlap_index;

	bool isReInitialize;
	bool isConvergence;
};

struct Species_Info
{
	int beginIndex;
	int length;
	double species_seed_fitness;
	int *index;
};
struct NewType
{
	int id;
	double data;
};
//bool Compare_NewType(NewType data1, NewType data2);
//bool Compare_NewType2(NewType data1, NewType data2);
bool Compare_NewType(NewType data1, NewType data2);
bool Compare_NewType_Des(NewType data1, NewType data2);
// ascending
bool Compare_NewType_Asc(NewType data1, NewType data2);


//const int timesOfRun = 51;
extern int timesOfRun;

//extern int exhaustSpace;	// strategy 2: mutate on its own axis range
//extern int timesToVisitBSP;	// strategy 1: mutate to nearest space
//extern int Mutation_Axis_Range;

extern int Mutate_Probability;
extern int Mutate_BSP_FULL;
extern int Mutate_On_Axis_Range;
extern int Mutate_To_Nearest_Space;
extern int No_Leaf_Initialization;
extern int Max_No_Visit;
extern int Converge_Times;

extern int BetterCounts;
extern int MiddleCounts;
extern int WorseCounts;

extern double DIST_INNER_CLUSTER;
extern double DIST_INTER_CLUSTER;

const double epsilon_set[] = { 1e-1,1e-2,1e-3,1e-4,1e-5 };
const int epsilon_set_size = 5;




double Fitness(double *particle, int &Fes, CEC2013 *pFunc);

void Popupation_Fitness(double **population, int population_size, int &Fes, double *result, CEC2013 *pFunc);

double Distance(double *vector1, vector<double> vector2, int dim);
double Distance(vector<double> vector1, vector<double> vector2, int dim);
double Distance(double *vector1, double *vector2, int dim);

int Close_Particle(double *child, double **population, int population_size, int dim);

void Get_Seeds(double **population, double *population_results, int population_size, int dim, vector<double> &seed_fitness, double radius);

int How_Many_Global_Optima(vector<double> seed_fitness, double epsilon, CEC2013 *pFunc);

int Compute_Global_Optima_Found(double **population, double *population_results, int population_size, int dim, double epsilon, CEC2013 *pFunc);
void Evolve(double **population, double**speed, double*population_result, double**personal_best,
	double*personal_best_result, int *neighbor_best_index, double *UBound, double *LBound,
	int population_size, double weight, double phi, int dim, int &Fes, CEC2013 *pFunc);

void Generate_Child(double *child, double *child_speed, double *exemplar1, double *exemplar2, double *UBound, double *LBound, double weight, double phi, int dim);

void printPop(double **population, double *population_results, int population_size, int dim);
void printArchive(vector<vector<double> >archive, vector<double> fitnessArchive, int dim);
// V8 
// add m, the number of individuals in one group that is nonoverlapping
void Get_Neighbor_Best_m_LHC(double *personal_best_results, int *neighbor_best_index, int population_size, int m);

//void Evolve_NrNSDE_V2_1(double **population, double **population_new, double*results, double *individual_trial,
//	double *individual_best, double &individual_best_result, double *UBound, double *LBound,
//	int population_size, int neighbor_size, double F, double CR, int dim, int &Fes, CEC2013 *pFunc,
//	boost::mt19937 &generator, double **species, double *species_results,
//	std::vector<std::vector<double> > &archive, std::vector<double> &fitnessArchive, double &maxFitness,
//	double *mutated_Chromosome, vector<Selection_Strategy> &Leaf_Vector,
//	NonRevisiting_GA *cur_GA, int axis_range, const double &radius);
//
//
////for V8
//
//
//
//
//void VisitStrategy_RPSO_V5(NonRevisiting_GA *cur_GA, int axis_range, double **population, double*population_result, double **speed, double**personal_best,
//	double*personal_best_result, int *neighbor_best_index, double *UBound, double *LBound,
//	std::vector<std::vector<double> > &archive, std::vector<double> &fitnessArchive, double &maxFitness,
//	int population_size, int dim, int &Fes, CEC2013 *pFunc,
//	double *mutated_Chromosome, boost::mt19937 &generator, vector<Selection_Strategy> &Leaf_Vector, const double &radius, const double &mutation_ratio,
//	int m);
//void VisitStrategy_V5_2(double **population, double **population_new, double*population_result, double *individual_trial,
//	double *individual_best, double &individual_best_result, double *UBound, double *LBound,
//	int population_size, int m, double F, double CR, int dim, int &Fes, CEC2013 *pFunc,
//	boost::mt19937 &generator, double **species, double *species_results,
//	std::vector<std::vector<double> > &archive, std::vector<double> &fitnessArchive, double &maxFitness,
//	double *mutated_Chromosome, vector<Selection_Strategy> &Leaf_Vector,
//	NonRevisiting_GA *cur_GA, int axis_range, const double &radius);
//void VisitStrategy_V5_1(double **population, double **population_new, double*population_result, double *individual_trial,
//	double *individual_best, double &individual_best_result, double *UBound, double *LBound,
//	int population_size, int m, double F, double CR, int dim, int &Fes, CEC2013 *pFunc,
//	boost::mt19937 &generator, double **species, double *species_results,
//	std::vector<std::vector<double> > &archive, std::vector<double> &fitnessArchive, double &maxFitness,
//	double *mutated_Chromosome, vector<Selection_Strategy> &Leaf_Vector,
//	NonRevisiting_GA *cur_GA, int axis_range, const double &radius);
//void VisitStrategy_RPSO_V5_wo_Inter(NonRevisiting_GA *cur_GA, int axis_range, double **population, double*population_result, double **speed, double**personal_best,
//	double*personal_best_result, int *neighbor_best_index, double *UBound, double *LBound,
//	std::vector<std::vector<double> > &archive, std::vector<double> &fitnessArchive, double &maxFitness,
//	int population_size, int dim, int &Fes, CEC2013 *pFunc,
//	double *mutated_Chromosome, boost::mt19937 &generator, vector<Selection_Strategy> &Leaf_Vector, const double &radius, const double &mutation_ratio,
//	int m);
//void VisitStrategy_V5_wo_Inter(double **population, double **population_new, double*population_result, double *individual_trial,
//	double *individual_best, double &individual_best_result, double *UBound, double *LBound,
//	int population_size, int m, double F, double CR, int dim, int &Fes, CEC2013 *pFunc,
//	boost::mt19937 &generator, double **species, double *species_results,
//	std::vector<std::vector<double> > &archive, std::vector<double> &fitnessArchive, double &maxFitness,
//	double *mutated_Chromosome, vector<Selection_Strategy> &Leaf_Vector,
//	NonRevisiting_GA *cur_GA, int axis_range, const double &radius);
//void VisitStrategy_V6_wo_Inter(double **population, double **population_new, double*population_result, double *individual_trial,
//	double *individual_best, double &individual_best_result, double *UBound, double *LBound,
//	int population_size, int m, double F, double CR, int dim, int &Fes, CEC2013 *pFunc,
//	boost::mt19937 &generator, double **species, double *species_results,
//	std::vector<std::vector<double> > &archive, std::vector<double> &fitnessArchive, double &maxFitness,
//	double *mutated_Chromosome, vector<Selection_Strategy> &Leaf_Vector,
//	NonRevisiting_GA *cur_GA, int axis_range, const double &inner_radius, const double &inter_radius);
////Éú³ÉÖÐ¼ä¸öÌå
//void VisitStrategy_V6(double **population, double **population_new, double*population_result, double *individual_trial,
//	double *individual_best, double &individual_best_result, double *UBound, double *LBound,
//	int population_size, int m, double F, double CR, int dim, int &Fes, CEC2013 *pFunc,
//	boost::mt19937 &generator, double **species, double *species_results,
//	std::vector<std::vector<double> > &archive, std::vector<double> &fitnessArchive, double &maxFitness,
//	double *mutated_Chromosome, vector<Selection_Strategy> &Leaf_Vector,
//	NonRevisiting_GA *cur_GA, int axis_range, const double &inner_radius, const double &inter_radius);
////¸ßË¹Éú³É
//void VisitStrategy_V6_1(double **population, double **population_new, double*population_result, double *individual_trial,
//	double *individual_best, double &individual_best_result, double *UBound, double *LBound,
//	int population_size, int m, double F, double CR, int dim, int &Fes, CEC2013 *pFunc,
//	boost::mt19937 &generator, double **species, double *species_results,
//	std::vector<std::vector<double> > &archive, std::vector<double> &fitnessArchive, double &maxFitness,
//	double *mutated_Chromosome, vector<Selection_Strategy> &Leaf_Vector,
//	NonRevisiting_GA *cur_GA, int axis_range, const double &inner_radius, const double &inter_radius);
////¸ßË¹Éú³É£¬¸ÅÂÊ½ÓÊÜ
//void VisitStrategy_V6_2(double **population, double **population_new, double*population_result, double *individual_trial,
//	double *individual_best, double &individual_best_result, double *UBound, double *LBound,
//	int population_size, int m, double F, double CR, int dim, int &Fes, CEC2013 *pFunc,
//	boost::mt19937 &generator, double **species, double *species_results,
//	std::vector<std::vector<double> > &archive, std::vector<double> &fitnessArchive, double &maxFitness,
//	double *mutated_Chromosome, vector<Selection_Strategy> &Leaf_Vector,
//	NonRevisiting_GA *cur_GA, int axis_range, const double &inner_radius, const double &inter_radius);
//// ×î¼Ñ¸öÌåµÄÖÐµã
//void VisitStrategy_V6_3(double **population, double **population_new, double*population_result, double *individual_trial,
//	double *individual_best, double &individual_best_result, double *UBound, double *LBound,
//	int population_size, int m, double F, double CR, int dim, int &Fes, CEC2013 *pFunc,
//	boost::mt19937 &generator, double **species, double *species_results,
//	std::vector<std::vector<double> > &archive, std::vector<double> &fitnessArchive, double &maxFitness,
//	double *mutated_Chromosome, vector<Selection_Strategy> &Leaf_Vector,
//	NonRevisiting_GA *cur_GA, int axis_range, const double &inner_radius, const double &inter_radius);
//void VisitStrategy_V6_4(double **population, double **population_new, double*population_result, double *individual_trial,
//	double *individual_best, double &individual_best_result, double *UBound, double *LBound,
//	int population_size, int m, double F, double CR, int dim, int &Fes, CEC2013 *pFunc,
//	boost::mt19937 &generator, double **species, double *species_results,
//	std::vector<std::vector<double> > &archive, std::vector<double> &fitnessArchive, double &maxFitness,
//	double *mutated_Chromosome, vector<Selection_Strategy> &Leaf_Vector,
//	NonRevisiting_GA *cur_GA, int axis_range, const double &inner_radius, const double &inter_radius);
//
//
//void VisitStrategy_V5(double **population, double **population_new, double*population_result, double *individual_trial,
//	double *individual_best, double &individual_best_result, double *UBound, double *LBound,
//	int population_size, int neighbor_size, double F, double CR, int dim, int &Fes, CEC2013 *pFunc,
//	boost::mt19937 &generator, double **species, double *species_results,
//	std::vector<std::vector<double> > &archive, std::vector<double> &fitnessArchive, double &maxFitness,
//	double *mutated_Chromosome, vector<Selection_Strategy> &Leaf_Vector,
//	NonRevisiting_GA *cur_GA, int axis_range, const double &radius);
//
//void Evolve_NrNSDE_V5_SDE(double **population, double **population_new, double*results, double *individual_trial,
//	double *individual_best, double &individual_best_result, double *UBound, double *LBound,
//	int population_size, int neighbor_size, double F, double CR, int dim, int &Fes, CEC2013 *pFunc,
//	boost::mt19937 &generator, double **species, double *species_results,
//	std::vector<std::vector<double> > &archive, std::vector<double> &fitnessArchive, double &maxFitness,
//	double *mutated_Chromosome, vector<Selection_Strategy> &Leaf_Vector,
//	NonRevisiting_GA *cur_GA, int axis_range, const double &radius);
//
////V7: ¶à¸ö¸öÌåµ½ÖÐÐÄµÄÅ·Ê½¾àÀë£¬ÇóÆ½¾ù
//void VisitStrategy_V6(double **population, double **population_new, double*population_result, double *individual_trial,
//	double *individual_best, double &individual_best_result, double *UBound, double *LBound,
//	int population_size, int m, double F, double CR, int dim, int &Fes, CEC2013 *pFunc,
//	boost::mt19937 &generator, double **species, double *species_results,
//	std::vector<std::vector<double> > &archive, std::vector<double> &fitnessArchive, double &maxFitness,
//	double *mutated_Chromosome, vector<Selection_Strategy> &Leaf_Vector,
//	NonRevisiting_GA *cur_GA, int axis_range, const double &inner_radius, const double &inter_radius);


void Evolve_NrNSDE_V5_random_centroid(double **population, double **population_new, double*results, double *individual_trial,
	double *individual_best, double &individual_best_result, double *UBound, double *LBound,
	int population_size, int neighbor_size, double F, double CR, int dim, int &Fes, CEC2013 *pFunc,
	boost::mt19937 &generator, double **species, double *species_results,
	std::vector<std::vector<double> > &archive, std::vector<double> &fitnessArchive, double &maxFitness,
	double *mutated_Chromosome, vector<Selection_Strategy> &Leaf_Vector,
	NonRevisiting_GA *cur_GA, int axis_range, const double &radius);
void Evolve_RmPSO_LHC(NonRevisiting_GA *cur_GA, int axis_range, double **population, double*population_result, double **speed, double**personal_best,
	double*personal_best_result, int *neighbor_best_index, double *UBound, double *LBound,
	std::vector<std::vector<double> > &archive, std::vector<double> &fitnessArchive, double &maxFitness,
	int population_size, int dim, int &Fes, CEC2013 *pFunc,
	double *mutated_Chromosome, boost::mt19937 &generator, vector<Selection_Strategy> &Leaf_Vector, const double &radius, const double &mutation_ratio,
	int m);

void VisitStrategy_V9(double **population, double **population_new, double*population_result, double *individual_trial,
	double *individual_best, double &individual_best_result, double *UBound, double *LBound,
	int population_size, int m, double F, double CR, int dim, int &Fes, CEC2013 *pFunc,
	boost::mt19937 &generator, double **species, double *species_results,
	std::vector<std::vector<double> > &archive, std::vector<double> &fitnessArchive, double &maxFitness,
	double *mutated_Chromosome, vector<Selection_Strategy> &Leaf_Vector,
	NonRevisiting_GA *cur_GA, int axis_range, const double &inner_radius, const double &inter_radius);
void VisitStrategy_V8_random_reInti(double **population, double **population_new, double*population_result, double *individual_trial,
	double *individual_best, double &individual_best_result, double *UBound, double *LBound,
	int population_size, int m, double F, double CR, int dim, int &Fes, CEC2013 *pFunc,
	boost::mt19937 &generator, double **species, double *species_results,
	std::vector<std::vector<double> > &archive, std::vector<double> &fitnessArchive, double &maxFitness,
	double *mutated_Chromosome, vector<Selection_Strategy> &Leaf_Vector,
	NonRevisiting_GA *cur_GA, int axis_range, const double &inner_radius, const double &inter_radius);
void VisitStrategy_V8_trace(double **population, double **population_new, double*population_result, double *individual_trial,
	double *individual_best, double &individual_best_result, double *UBound, double *LBound,
	int population_size, int m, double F, double CR, int dim, int &Fes, CEC2013 *pFunc,
	boost::mt19937 &generator, double **species, double *species_results,
	std::vector<std::vector<double> > &archive, std::vector<double> &fitnessArchive, double &maxFitness,
	double *mutated_Chromosome, vector<Selection_Strategy> &Leaf_Vector,
	NonRevisiting_GA *cur_GA, int axis_range, const double &inner_radius, const double &inter_radius,
	int func, int *traceFEs, int &trace_pos, const int &times);


void VisitStrategy_V8_RPSO(double **population, double*population_result,
	double **speed, double**personal_best, double*personal_best_result, int *neighbor_best_index,
	double *UBound, double *LBound, int population_size, int m, int dim, int &Fes, CEC2013 *pFunc,
	boost::mt19937 &generator, std::vector<std::vector<double> > &archive, std::vector<double> &fitnessArchive, double &maxFitness,
	double *mutated_Chromosome, vector<Selection_Strategy> &Leaf_Vector,
	NonRevisiting_GA *cur_GA, int axis_range, const double &inner_radius, const double &inter_radius);
void VisitStrategy_V8_RPSO_random_reInti(double **population, double*population_result,
	double **speed, double**personal_best, double*personal_best_result, int *neighbor_best_index,
	double *UBound, double *LBound, int population_size, int m, int dim, int &Fes, CEC2013 *pFunc,
	boost::mt19937 &generator, std::vector<std::vector<double> > &archive, std::vector<double> &fitnessArchive, double &maxFitness,
	double *mutated_Chromosome, vector<Selection_Strategy> &Leaf_Vector,
	NonRevisiting_GA *cur_GA, int axis_range, const double &inner_radius, const double &inter_radius);
#endif // SELF_DEFINE_FUNCTIONS_H_INCLUDED