/*
 * VisitStrategy.h
 *
 *  Created on: Aug 2, 2017
 *      Author: hting
 */

#ifndef VISITSTRATEGY_H_
#define VISITSTRATEGY_H_

#include "../../header.h"
#include "../../Archive/Archive.h"
#include "../../BSP/ACC_Math_Kernel/ACC_Math_Kernel.h"
#include "../../BSP/GeneticAlgorithm_Kernel/GeneticAlgorithm_Kernel.h"
#include "../../BSP/NonRevisitingGA_Kernel/NonRevisitingGA_Kernel.h"
#include "../../BSP/Optimization_TestFunction_Kernel/Optimization_TestFunction_Kernel.h"
#include "../../BSP/RecordRevisitingScheme_Kernel/RecordRevisitingScheme_Kernel.h"
#include "../../BSP/SelectionStrategey_Kernel/Selection_Strategy.h"

int VisitStrateg_Record(double *cur_Chromosome, double cur_Chromonsome_Fitness, int index,
	NonRevisiting_GA *cur_GA, int *Current_ID, double **cur_Interval, vector<Selection_Strategy> &Select_Vector);
double Cal_Mean_Dist(double **population, int begin_index, int dim, int length);
int VisitStrategy_Reintialization(double *cur_Chromosome, double *mutated_Chromosome,
	int visitTimes, int index, NonRevisiting_GA *cur_GA, int axis_range,
	int dim, double *LBound, double *UBound, boost::mt19937 &generator, vector<Selection_Strategy> &Leaf_Vector, const double &mutation_ratio,
	const int Fes, const int maxFes);
void R3PSO_LHC_VisitStrategy(NonRevisiting_GA *cur_GA, int axis_range, double **population, double*population_result, double **speed, double**personal_best,
	double*personal_best_result, int *neighbor_best_index, double *UBound, double *LBound,
	std::vector<std::vector<double> > &archive, std::vector<double> &fitnessArchive, double &maxFitness,
	int population_size, int dim, int &Fes, CEC2013 *pFunc,
	double *mutated_Chromosome, boost::mt19937 &generator, vector<Selection_Strategy> &Leaf_Vector, const double &radius, const double &mutation_ratio);
void R2PSO_LHC_VisitStrategy(NonRevisiting_GA *cur_GA, int axis_range, double **population, double*population_result, double **speed, double**personal_best,
	double*personal_best_result, int *neighbor_best_index, double *UBound, double *LBound,
	std::vector<std::vector<double> > &archive, std::vector<double> &fitnessArchive, double &maxFitness,
	int population_size, int dim, int &Fes, CEC2013 *pFunc,
	double *mutated_Chromosome, boost::mt19937 &generator, vector<Selection_Strategy> &Leaf_Vector, const double &radius, const double &mutation_ratio);


void R3PSO_LHC_VisitStrategy_Random(NonRevisiting_GA *cur_GA, int axis_range, double **population, double*population_result, double **speed, double**personal_best,
	double*personal_best_result, int *neighbor_best_index, double *UBound, double *LBound,
	std::vector<std::vector<double> > &archive, std::vector<double> &fitnessArchive, double &maxFitness,
	int population_size, int dim, int &Fes, CEC2013 *pFunc,
	double *mutated_Chromosome, boost::mt19937 &generator, vector<Selection_Strategy> &Leaf_Vector, const double &radius, const double &mutation_ratio);
void R2PSO_LHC_VisitStrategy_Random(NonRevisiting_GA *cur_GA, int axis_range, double **population, double*population_result, double **speed, double**personal_best,
	double*personal_best_result, int *neighbor_best_index, double *UBound, double *LBound,
	std::vector<std::vector<double> > &archive, std::vector<double> &fitnessArchive, double &maxFitness,
	int population_size, int dim, int &Fes, CEC2013 *pFunc,
	double *mutated_Chromosome, boost::mt19937 &generator, vector<Selection_Strategy> &Leaf_Vector, const double &radius, const double &mutation_ratio);


#endif /* VISITSTRATEGY_H_ */
