/******************************************************************************************/
//	Program Name	:	Non-Revisiting Genetic Algorithm (NrGA)
//	File Name		:	Optimization_TestFunction_Kernel.c
//	Auther			:	Dr. Chow Chi Kin, Dr. Yuen Shiu Yin
//	Edit by			:	Leung Shing Wa
//	University		:	City University of Hong Kong
//	Department		:	Electronic Engineering
//	Last Update		:	10 Sep 2009
//	Reference		:	A Genetic Algorithm that Adaptively Mutates and Never Revisits, 
//						IEEE Transactions on Evolutionary Computation, 
//						Vol 13(2) (April 2009) 454-472.
//
//	Discription		:	All Optimization Test functions can be found in here.

#include "../../BSP/Optimization_TestFunction_Kernel/Optimization_TestFunction_Kernel.h"

#include<stdio.h>
#include<stdlib.h>
#include<malloc.h>
#include<math.h>
#include<time.h>

#include "../../BSP/ACC_Math_Kernel/ACC_Math_Kernel.h"

// Global variables
int Optimization_TestFunction_RotationFlag;
int Optimization_TestFunction_Dimension;
double Optimization_TestFunction_AxisRange;

int Optimization_TestFunction_RealRoyalRoadTau;
int Optimization_TestFunction_MultiRealValue_BlockSize;

double Optimization_TestFunction_LinearWeigth[256];
double Optimization_TestFunction_XORTReeNode[2][256];
double *Optimization_TestFunction_TruthTableValue;

double Optimization_TestFunction_Unbounded_Bias;

double Optimization_TestFunction_NoiseLevel;
double (*Optimization_TestFunction_TruthTableFunction)(double*);
double (*Optimization_TestFunction_Pointer)(double*);

// original file

double **Optimization_TestFunction_RotationMatrix, *Scaled_Solution, *Rotated_Solution;

double Optimization_TestFunction_Foxholes_a1[] = {-32.0, -16.0, 0, 16.0, 32.0, -32.0,
										-32.0, -16.0, 0, 16.0, 32.0, -32.0,
										-32.0, -16.0, 0, 16.0, 32.0, -32.0,
										-32.0, -16.0, 0, 16.0, 32.0, -32.0,
										-32.0, -16.0, 0, 16.0, 32.0, -32.0};
double Optimization_TestFunction_Foxholes_a2[] = {-32.0, -32.0, -32.0, -32.0, -32.0,
										-16.0, -16.0, -16.0, -16.0, -16.0,
										0.0, 0.0, 0.0, 0.0, 0.0,
										16.0, 16.0, 16.0, 16.0, 16.0,
										32.0, 32.0, 32.0, 32.0, 32.0};

double Optimization_TestFunction_OddSquare_b[] = {1.0, 1.3, 0.8, -0.4, -1.3,
												  1.6, -2.0, -6.0, 0.5, 1.4};

void Optimization_TestFunction_Manager_RealValued(int Optimization_TestFunction_Index, char *Optimization_TestFunction_Name)
{
	switch(Optimization_TestFunction_Index)
	{
		case 0:			//------------ Sphereical
			Optimization_TestFunction_Pointer = Optimization_TestFunction_Sphereical;
			sprintf(Optimization_TestFunction_Name, "Spherical");
			break;
		case 1:			//------------ Schwefel 2.22
			Optimization_TestFunction_Pointer = Optimization_TestFunction_Schwefel2_22;
			sprintf(Optimization_TestFunction_Name, "Schwefel2-22");
			break;
		case 2:			//------------ Schwefel 1.02
			Optimization_TestFunction_Pointer = Optimization_TestFunction_Schwefel1_02;
			sprintf(Optimization_TestFunction_Name, "Schwefel1-02");
			break;
		case 3:			//------------ Schwefel 2.21
			Optimization_TestFunction_Pointer = Optimization_TestFunction_Schwefel2_21;
			sprintf(Optimization_TestFunction_Name, "Schwefel2-21");
			break;
		case 4:			//------------ Rosenbrock
			Optimization_TestFunction_Pointer = Optimization_TestFunction_Rosenbrock;
			sprintf(Optimization_TestFunction_Name, "Rosenbrock");
			break;
		case 5:			//------------ Quartic
			Optimization_TestFunction_Pointer = Optimization_TestFunction_Quartic;
			sprintf(Optimization_TestFunction_Name, "Quartic");
			break;
		case 6:			//------------ Generalized Rastrigin
			Optimization_TestFunction_Pointer = Optimization_TestFunction_GeneralizedRastrigin;
			sprintf(Optimization_TestFunction_Name, "GeneralizedRastrigin");
			break;
		case 7:			//------------ Generalized Griewank
			Optimization_TestFunction_Pointer = Optimization_TestFunction_GeneralizedGriewank;
			sprintf(Optimization_TestFunction_Name, "GeneralizedGriewank");
			break;
		case 8:			//------------ Schwefel 2.26
			Optimization_TestFunction_Pointer = Optimization_TestFunction_Schwefel2_26;
			sprintf(Optimization_TestFunction_Name, "Schwefel2-26");
			break;
		case 9:			//------------ Ackley
			Optimization_TestFunction_Pointer = Optimization_TestFunction_Ackley;
			sprintf(Optimization_TestFunction_Name, "Ackley");
			break;
		case 10:		//------------ Foxholes
			Optimization_TestFunction_Pointer = Optimization_TestFunction_Foxholes;
			sprintf(Optimization_TestFunction_Name, "Foxholes");
			break;
		case 11:		//------------ Six-Hump Camel-Back
			Optimization_TestFunction_Pointer = Optimization_TestFunction_SixHumpCamelBack;
			sprintf(Optimization_TestFunction_Name, "SixHumpCamelBack");
			break;
		case 12:		//------------ Branin
			Optimization_TestFunction_Pointer = Optimization_TestFunction_Branin;
			sprintf(Optimization_TestFunction_Name, "Branin");
			break;
		case 13:		//------------ Goldstein-Price
			Optimization_TestFunction_Pointer = Optimization_TestFunction_GoldsteinPrice;
			sprintf(Optimization_TestFunction_Name, "GoldsteinPrice");
			break;
		case 14:		//------------ High Conditioned Elliptic
			Optimization_TestFunction_Pointer = Optimization_TestFunction_HighConditionedElliptic;
			sprintf(Optimization_TestFunction_Name, "HighConditionedElliptic");
			break;
		case 15:		//------------ Weierstrass
			Optimization_TestFunction_Pointer = Optimization_TestFunction_Weierstrass;
			sprintf(Optimization_TestFunction_Name, "Weierstrass");
			break;
		case 16:		//------------ Hybrid Composition 1
			Optimization_TestFunction_Pointer = Optimization_TestFunction_HybridComposition01;
			sprintf(Optimization_TestFunction_Name, "HybridComposition01");
			break;
		case 17:		//------------ Levy
			Optimization_TestFunction_Pointer = Optimization_TestFunction_Levy;
			sprintf(Optimization_TestFunction_Name, "Levy");
			break;
		case 18:		//------------ Pern
			Optimization_TestFunction_Pointer = Optimization_TestFunction_Pern;
			sprintf(Optimization_TestFunction_Name, "Pern");
			break;
		case 19:		//------------ Zakharov
			Optimization_TestFunction_Pointer = Optimization_TestFunction_Zakharov;
			sprintf(Optimization_TestFunction_Name, "Zakharov");
			break;
		case 20:		//------------ Alpine
			Optimization_TestFunction_Pointer = Optimization_TestFunction_Alpine;
			sprintf(Optimization_TestFunction_Name, "Alpine");
			break;
		case 21:		//------------ Pathological
			Optimization_TestFunction_Pointer = Optimization_TestFunction_Pathological;
			sprintf(Optimization_TestFunction_Name, "Pathological");
			break;
		case 22:		//------------ InvertedCosineWave
			Optimization_TestFunction_Pointer = Optimization_TestFunction_InvertedCosineMixture;
			sprintf(Optimization_TestFunction_Name, "InvertedCosineMixture");
			break;
		case 23:		//------------ InvertedCosineMixture
			Optimization_TestFunction_Pointer = Optimization_TestFunction_InvertedCosineWave;
			sprintf(Optimization_TestFunction_Name, "InvertedCosineWave");
			break;
		case 24:		//------------ EpistaticMichalewicz
			Optimization_TestFunction_Pointer = Optimization_TestFunction_EpistaticMichalewicz;
			sprintf(Optimization_TestFunction_Name, "EpistaticMichalewicz");
			break;
		case 25:		//------------ LevyMontalvo
			Optimization_TestFunction_Pointer = Optimization_TestFunction_LevyMontalvo;
			sprintf(Optimization_TestFunction_Name, "LevyMontalvo");
			break;
		case 26:		//------------ Neumaier3
			Optimization_TestFunction_Pointer = Optimization_TestFunction_Neumaier3;
			sprintf(Optimization_TestFunction_Name, "Neumaier3");
			break;
		case 27:		//------------ OddSquare
			Optimization_TestFunction_Pointer = Optimization_TestFunction_OddSquare;
			sprintf(Optimization_TestFunction_Name, "OddSquare");
			break;
		case 28:		//------------ Paviani
			Optimization_TestFunction_Pointer = Optimization_TestFunction_Paviani;
			sprintf(Optimization_TestFunction_Name, "Paviani");
			break;
		case 29:		//------------ Periodic
			Optimization_TestFunction_Pointer = Optimization_TestFunction_Periodic;
			sprintf(Optimization_TestFunction_Name, "Periodic");
			break;
		case 30:		//------------ Salomon
			Optimization_TestFunction_Pointer = Optimization_TestFunction_Salomon;
			sprintf(Optimization_TestFunction_Name, "Salomon");
			break;
		case 31:		//------------ Shubert
			Optimization_TestFunction_Pointer = Optimization_TestFunction_Shubert;
			sprintf(Optimization_TestFunction_Name, "Shubert");
			break;
		case 32:		//------------ Sinusoidal
			Optimization_TestFunction_Pointer = Optimization_TestFunction_Sinusoidal;
			sprintf(Optimization_TestFunction_Name, "Sinusoidal");
			break;
		case 33:		//------------ Michalewicz
			Optimization_TestFunction_Pointer = Optimization_TestFunction_Michalewicz;
			sprintf(Optimization_TestFunction_Name, "Michalewicz");
			break;
		case 34:		//------------ Whitely
			Optimization_TestFunction_Pointer = Optimization_TestFunction_Whitely;
			sprintf(Optimization_TestFunction_Name, "Whitely");
			break;
		case 35:		//------------ Easom
			Optimization_TestFunction_Pointer = Optimization_TestFunction_Easom;
			sprintf(Optimization_TestFunction_Name, "Easom");
			break;
	}

	return;
}

void Optimization_TestFunction_Manager_PesudoBoolean(int Optimization_TestFunction_Index, char *Optimization_TestFunction_Name)
{
	switch(Optimization_TestFunction_Index)
	{
		case 0:		//------------ Max-One
			Optimization_TestFunction_Pointer = Optimization_TestFunction_MaxOne;
			sprintf(Optimization_TestFunction_Name, "MaxOne");
			break;
		case 1:		//------------ Real Royal Road
			Optimization_TestFunction_Pointer = Optimization_TestFunction_RealRoyalRoad;
//			GA_TESTFUNCTION = Optimization_TestFunction_InverseRealRoyalRoad2001;
			sprintf(Optimization_TestFunction_Name, "RealRoyalRoad");
			Optimization_TestFunction_RealRoyalRoadTau = (int) (Optimization_TestFunction_Dimension * 0.1);
			break;
		case 2:		//------------ Linear
			Optimization_TestFunction_Pointer = Optimization_TestFunction_Linear;
			sprintf(Optimization_TestFunction_Name, "Linear");
			break;
		case 3:		//------------ Almost Positive
			Optimization_TestFunction_Pointer = Optimization_TestFunction_AlmostPositive;
			sprintf(Optimization_TestFunction_Name, "AlmostPositive");
			break;
		case 4:		//------------ Inverse Pseudo
			Optimization_TestFunction_Pointer = Optimization_TestFunction_InversePseudoModular;
			sprintf(Optimization_TestFunction_Name, "InversePseudo");
			break;
		case 5:		//------------ Short Path Constant
			Optimization_TestFunction_Pointer = Optimization_TestFunction_ShortPathConstant;
			sprintf(Optimization_TestFunction_Name, "ShortPathConstant");
			break;
		case 6:		//------------ Real Value
			Optimization_TestFunction_Pointer = Optimization_TestFunction_RealValue;
			sprintf(Optimization_TestFunction_Name, "RealValue");
			break;
		case 7:		//------------ XOR-Tree
			Optimization_TestFunction_Pointer = Optimization_TestFunction_XORTree;
			sprintf(Optimization_TestFunction_Name, "XORTree");
			break;
	}

	return;
}

void Optimization_TestFunction_Information(double **min_X, double **dX, double **TrueOptima, double *TrueOptimalFitness, int testfunction_index)
{
	int i;

	switch(testfunction_index)
	{
		case 0:
			for(i=0; i< Optimization_TestFunction_Dimension; ++i)
			{
				(*min_X)[i] = -10.0;
				(*dX)[i] = 20.0;
				(*TrueOptima)[i] = 0.0;
			}
			(*TrueOptimalFitness) = 0.0;
			break;
		case 1:
			for(i=0; i< Optimization_TestFunction_Dimension; ++i)
			{
				(*min_X)[i] = -10.0;
				(*dX)[i] = 20.0;
				(*TrueOptima)[i] = 0.0;
			}
			(*TrueOptimalFitness) = 0.0;
			break;
		case 2:
			for(i=0; i< Optimization_TestFunction_Dimension; ++i)
			{
				(*min_X)[i] = -100.0;
				(*dX)[i] = 200.0;
				(*TrueOptima)[i] = 0.0;
			}
			(*TrueOptimalFitness) = 0.0;
			break;
		case 3:
			for(i=0; i< Optimization_TestFunction_Dimension; ++i)
			{
				(*min_X)[i] = -29.0;
				(*dX)[i] = 60.0;
				(*TrueOptima)[i] = 0.0;
			}
			(*TrueOptimalFitness) = 0.0;
			break;
		case 4:
			for(i=0; i< Optimization_TestFunction_Dimension; ++i)
			{
				(*min_X)[i] = -1.28;
				(*dX)[i] = 2.56;
				(*TrueOptima)[i] = 1.0;
			}
			(*TrueOptimalFitness) = 0.0;
			break;
		case 5:
			for(i=0; i< Optimization_TestFunction_Dimension; ++i)
			{
				(*min_X)[i] = -5.12;
				(*dX)[i] = 10.24;
				(*TrueOptima)[i] = 0.0;
			}
			(*TrueOptimalFitness) = 0.0;
			break;
		case 6:
			for(i=0; i< Optimization_TestFunction_Dimension; ++i)
			{
				(*min_X)[i] = -600.0;
				(*dX)[i] = 1200.0;
				(*TrueOptima)[i] = 0.0;
			}
			(*TrueOptimalFitness) = 0.0;
			break;
		case 7:
			for(i=0; i< Optimization_TestFunction_Dimension; ++i)
			{
				(*min_X)[i] = -500.0;
				(*dX)[i] = 1000.0;
				(*TrueOptima)[i] = 0.0;
			}
			(*TrueOptimalFitness) = 0.0;
			break;
		case 8:
			for(i=0; i< Optimization_TestFunction_Dimension; ++i)
			{
				(*min_X)[i] = -32.0;
				(*dX)[i] = 64.0;
				(*TrueOptima)[i] = 420.9687;
			}
			(*TrueOptimalFitness) = -418.9829 * Optimization_TestFunction_Dimension;
			break;
		case 9:
			for(i=0; i< Optimization_TestFunction_Dimension; ++i)
			{
				(*min_X)[i] = -98.0;
				(*dX)[i] = 132.0;
				(*TrueOptima)[i] = 0.0;
			}
			(*TrueOptimalFitness) = 0.0;
			break;
		case 10:
			(*min_X)[0] = -98;
			(*min_X)[1] = -98;
			(*dX)[0] = 132.0;
			(*dX)[1] = 132.0;
			(*TrueOptima)[0] = -32.0;
			(*TrueOptima)[1] = -32.0;
			(*TrueOptimalFitness) = INF;
			break;
		case 11:
			(*min_X)[0] = -4.91017;
			(*min_X)[1] = -5.7126;
			(*dX)[0] = 10.0;
			(*dX)[1] = 10.0;
			(*TrueOptima)[0] = 0.0893;
			(*TrueOptima)[1] = -0.7126;
			(*TrueOptimalFitness) = INF;
			break;
		case 12:
			(*min_X)[0] = -8.142;
			(*min_X)[1] = -12.275;
			(*dX)[0] = 15.0;
			(*dX)[1] = 15.0;
			(*TrueOptima)[0] = -3.142;
			(*TrueOptima)[1] = 12.275;
			(*TrueOptimalFitness) = INF;
			break;
		case 13:
			(*min_X)[0] = -2;
			(*min_X)[1] = -2;
			(*dX)[0] = 4.0;
			(*dX)[1] = 4.0;
			(*TrueOptima)[0] = 0.0;
			(*TrueOptima)[1] = -1.0;
			(*TrueOptimalFitness) = INF;
			break;
		case 14:      //============ High Conditioned Elliptic Function
			for(i=0; i< Optimization_TestFunction_Dimension; ++i)
			{
				(*min_X)[i] = -100.0;
				(*dX)[i] = 200.0;
				(*TrueOptima)[i] = -100.0;
			}
			(*TrueOptimalFitness) = 0.0;
			break;
		case 15:	//============ Weierstrass Function
			for(i=0; i< Optimization_TestFunction_Dimension; ++i)
			{
				(*min_X)[i] = -0.5;
				(*dX)[i] = 1.0;
				(*TrueOptima)[i] = 0.0;
			}
			(*TrueOptimalFitness) = 0.0;
			break;
		case 16:	//============ Hybrid Composition 1
			for(i=0; i< Optimization_TestFunction_Dimension; ++i)
			{
				(*min_X)[i] = -5.0;
				(*dX)[i] = 10.0;
				(*TrueOptima)[i] = 0.0;
			}
			(*TrueOptimalFitness) = 4500.0;
			break;
		case 17:	//============ Levy
			for(i=0; i< Optimization_TestFunction_Dimension; ++i)
			{
				(*min_X)[i] = -10.0;
				(*dX)[i] = 20.0;
				(*TrueOptima)[i] = 1.0;
			}
			(*TrueOptimalFitness) = 0.0;
			break;
		case 18:	//============ Pern
			for(i=0; i< Optimization_TestFunction_Dimension; ++i)
			{
				(*min_X)[i] = -Optimization_TestFunction_Dimension;
				(*dX)[i] = 2.0 * Optimization_TestFunction_Dimension;
				(*TrueOptima)[i] = i+1;
			}
			(*TrueOptimalFitness) = INF;
			break;
		case 19:	//============ Zakharov
			for(i=0; i< Optimization_TestFunction_Dimension; ++i)
			{
				(*min_X)[i] = -5.0;
				(*dX)[i] = 15.0;
				(*TrueOptima)[i] = 0.0;
			}
			(*TrueOptimalFitness) = 0.0;
			break;
		case 20:	//============ Alpine
			for(i=0; i< Optimization_TestFunction_Dimension; ++i)
			{
				(*min_X)[i] = -10.0;
				(*dX)[i] = 20.0;
				(*TrueOptima)[i] = 0.0;
			}
			(*TrueOptimalFitness) = 0.0;
			break;
		case 21:	//============ Pathological
			for(i=0; i< Optimization_TestFunction_Dimension; ++i)
			{
				(*min_X)[i] = -100.0;
				(*dX)[i] = 200.0;
				(*TrueOptima)[i] = 0.0;
			}
			(*TrueOptimalFitness) = 0.0;
			break;
		case 22:	//============ InvertedCosineWave
			for(i=0; i< Optimization_TestFunction_Dimension; ++i)
			{
				(*min_X)[i] = -5.0;
				(*dX)[i] = 10.0;
				(*TrueOptima)[i] = 0.0;
			}
			(*TrueOptimalFitness) = -Optimization_TestFunction_Dimension + 1.0;
			break;
		case 23:	//============ CosineMixture
			for(i=0; i< Optimization_TestFunction_Dimension; ++i)
			{
				(*min_X)[i] = -1.0;
				(*dX)[i] = 2.0;
				(*TrueOptima)[i] = 0.0;
			}
			(*TrueOptimalFitness) = 0.0;
			break;
		case 24:	//============ EpistaticMichalewicz
			for(i=0; i< Optimization_TestFunction_Dimension; ++i)
			{
				(*min_X)[i] = 0.0;
				(*dX)[i] = PI;
				(*TrueOptima)[i] = 0.0;
			}
			(*TrueOptimalFitness) = INF;
			break;
		case 25:	//============ LevyMontalvo
			for(i=0; i< Optimization_TestFunction_Dimension; ++i)
			{
				(*min_X)[i] = -5.0;
				(*dX)[i] = 10.0;
				(*TrueOptima)[i] = 1.0;
			}
			(*TrueOptimalFitness) = 0.0;
			break;
		case 26:	//============ Neumaier3
			for(i=0; i< Optimization_TestFunction_Dimension; ++i)
			{
				(*min_X)[i] = -pow(Optimization_TestFunction_Dimension, 2);
				(*dX)[i] = 2.0 * pow(Optimization_TestFunction_Dimension, 2);
				(*TrueOptima)[i] = (i+1.0) * (Optimization_TestFunction_Dimension + i);
			}
			(*TrueOptimalFitness) = -Optimization_TestFunction_Dimension * (Optimization_TestFunction_Dimension+4.0) * (Optimization_TestFunction_Dimension-1.0) / 6.0;
			break;
		case 27:	//============ OddSquare
			for(i=0; i< Optimization_TestFunction_Dimension; ++i)
			{
				(*min_X)[i] = -15.0;
				(*dX)[i] = 30.0;
				(*TrueOptima)[i] = 0.0;
			}
			(*TrueOptimalFitness) = INF;
			break;
		case 28:	//============ Paviani
			for(i=0; i< Optimization_TestFunction_Dimension; ++i)
			{
				(*min_X)[i] = 2.0;
				(*dX)[i] = 8.0;
				(*TrueOptima)[i] = 0.0;
			}
			(*TrueOptimalFitness) = INF;
			break;
		case 29:	//============ Periodic
			for(i=0; i< Optimization_TestFunction_Dimension; ++i)
			{
				(*min_X)[i] = -10.0;
				(*dX)[i] = 20.0;
				(*TrueOptima)[i] = 0.0;
			}
			(*TrueOptimalFitness) = 0.9;
			break;
		case 30:	//============ Salomon
			for(i=0; i< Optimization_TestFunction_Dimension; ++i)
			{
				(*min_X)[i] = -100.0;
				(*dX)[i] = 200.0;
				(*TrueOptima)[i] = 0.0;
			}
			(*TrueOptimalFitness) = 0.0;
			break;
		case 31:	//============ Shubert
			for(i=0; i< Optimization_TestFunction_Dimension; ++i)
			{
				(*min_X)[i] = -10.0;
				(*dX)[i] = 20.0;
				(*TrueOptima)[i] = 0.0;
			}
			(*TrueOptimalFitness) = INF;
			break;
		case 32:	//============ Sinusoidal
			for(i=0; i< Optimization_TestFunction_Dimension; ++i)
			{
				(*min_X)[i] = -10.0;
				(*dX)[i] = 20.0;
				(*TrueOptima)[i] = 0.0;
			}
			(*TrueOptimalFitness) = INF;
			break;
		case 33:	//============ Michalewicz
			for(i=0; i< Optimization_TestFunction_Dimension; ++i)
			{
				(*min_X)[i] = 0.0;
				(*dX)[i] = 180.0;
				(*TrueOptima)[i] = 90.0 + 30.0;
			}
			(*TrueOptimalFitness) = -3.5;
			break;
		case 34:	//============ Whitely
			for(i=0; i< Optimization_TestFunction_Dimension; ++i)
			{
				(*min_X)[i] = -100.0;
				(*dX)[i] = 200.0;
				(*TrueOptima)[i] = 0.0;
			}
			(*TrueOptimalFitness) = INF;
			break;
		case 35:	//============ Whitely
			for(i=0; i< Optimization_TestFunction_Dimension; ++i)
			{
				(*min_X)[i] = -10.0;
				(*dX)[i] = 20.0;
				(*TrueOptima)[i] = PI;
			}
			(*TrueOptimalFitness) = -1.0;
			break;
	}

	return;
}

//================================================================================================================
//==================================== Real Value Objective Function =============================================
//================================================================================================================

double Optimization_TestFunction_Planar(double *solution)
{
	double fitness_value;

	Optimization_TestFunction_RotationSolution(solution, 1.0, 2.0);

	fitness_value = Rotated_Solution[0];

	if(Optimization_TestFunction_NoiseLevel >= 0.0 && Optimization_TestFunction_NoiseLevel <= 1000.0)
		fitness_value += randn(0, Optimization_TestFunction_NoiseLevel);

	return fitness_value;
}

double Optimization_TestFunction_Sphereical(double *solution)
{
	int i;
	double fitness_value;

	Optimization_TestFunction_RotationSolution(solution, 10.0, 20.0);

	fitness_value = 0.0;
	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
		fitness_value += pow(Rotated_Solution[i], 2);

	if(Optimization_TestFunction_NoiseLevel >= 0.0 && Optimization_TestFunction_NoiseLevel <= 1000.0)
		fitness_value += randn(0, Optimization_TestFunction_NoiseLevel);

	return fitness_value;
}

double Optimization_TestFunction_Schwefel2_22(double *solution)
{
	int i;
	double sum_x, prd_x, fitness_value;

	Optimization_TestFunction_RotationSolution(solution, 10.0, 20.0);

	sum_x = 0.0;
	prd_x = 1.0;

	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
	{
		sum_x += fabs(Rotated_Solution[i]);
		prd_x *= fabs(Rotated_Solution[i]);
	}

	fitness_value = sum_x + prd_x;

	if(Optimization_TestFunction_NoiseLevel >= 0.0 && Optimization_TestFunction_NoiseLevel <= 1000.0)
		fitness_value += randn(0, Optimization_TestFunction_NoiseLevel);

	return fitness_value;
}

double Optimization_TestFunction_Schwefel1_02(double *solution)
{
	int i, j;
	double sum_x, fitness_value;

	Optimization_TestFunction_RotationSolution(solution, 100.0, 200.0);

	fitness_value = 0.0;
	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
	{
		sum_x = 0.0;
		for(j=0; j <=i; ++j)
			sum_x += Rotated_Solution[j];
		fitness_value += pow(sum_x, 2);
	}

	if(Optimization_TestFunction_NoiseLevel >= 0.0 && Optimization_TestFunction_NoiseLevel <= 1000.0)
		fitness_value += randn(0, Optimization_TestFunction_NoiseLevel);

	return fitness_value;
}

double Optimization_TestFunction_Schwefel2_21(double *solution)
{
	int i;
	double fitness_value;

	Optimization_TestFunction_RotationSolution(solution, 100.0, 200.0);

	fitness_value = fabs(Rotated_Solution[0]);
	for(i=1; i< Optimization_TestFunction_Dimension; ++i)
	{
		if(fitness_value < fabs(Rotated_Solution[i]))
			fitness_value = fabs(Rotated_Solution[i]);
	}

	if(Optimization_TestFunction_NoiseLevel >= 0.0 && Optimization_TestFunction_NoiseLevel <= 1000.0)
		fitness_value += randn(0, Optimization_TestFunction_NoiseLevel);

	return fitness_value;
}

double Optimization_TestFunction_Rosenbrock(double *solution)
{
	int i;
	double x1, x2, fitness_value;

	Optimization_TestFunction_RotationSolution(solution, 29.0, 60.0);

	fitness_value = 0.0;
	for(i=0; i< Optimization_TestFunction_Dimension-1; ++i)
	{
		x1 = Rotated_Solution[i];
		x2 = Rotated_Solution[i+1];

		fitness_value += 100.0 * pow(x2 - pow(x1, 2), 2) + pow(x1 - 1, 2);
	}

	if(Optimization_TestFunction_NoiseLevel >= 0.0 && Optimization_TestFunction_NoiseLevel <= 1000.0)
		fitness_value += randn(0, Optimization_TestFunction_NoiseLevel);

	return fitness_value;
}

double Optimization_TestFunction_Quartic(double *solution)
{
	int i;
	double fitness_value;

	Optimization_TestFunction_RotationSolution(solution, 1.28, 2.56);

	fitness_value = 0.0;
	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
		fitness_value += (i + 1) * pow(Rotated_Solution[i], 4) + (double) (rand() % 10000) / 10000.0;

	if(Optimization_TestFunction_NoiseLevel >= 0.0 && Optimization_TestFunction_NoiseLevel <= 1000.0)
		fitness_value += randn(0, Optimization_TestFunction_NoiseLevel);

	return fitness_value;
}

double Optimization_TestFunction_GeneralizedRastrigin(double *solution)
{
	int i;
	double fitness_value;

	Optimization_TestFunction_RotationSolution(solution, 5.12, 10.24);

	fitness_value = 0.0;
	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
	{
		fitness_value += (pow(Rotated_Solution[i], 2) -
						  10.0 * cos(2.0 * PI * Rotated_Solution[i]) +
						  10.0);
	}

	if(Optimization_TestFunction_NoiseLevel >= 0.0 && Optimization_TestFunction_NoiseLevel <= 1000.0)
		fitness_value += randn(0, Optimization_TestFunction_NoiseLevel);

	return fitness_value;
}

double Optimization_TestFunction_GeneralizedGriewank(double *solution)
{
	int i;
	double sum_x, product_cosx;
	double fitness_value;

	Optimization_TestFunction_RotationSolution(solution, 600.0, 1200.0);

	sum_x = 0.0;
	product_cosx = 1.0;
	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
	{
		sum_x += pow(Rotated_Solution[i], 2);
		product_cosx *= cos(Rotated_Solution[i] / sqrt(i+1));
	}

	fitness_value = sum_x/4000.0 - product_cosx + 1.0;

	if(Optimization_TestFunction_NoiseLevel >= 0.0 && Optimization_TestFunction_NoiseLevel <= 1000.0)
		fitness_value += randn(0, Optimization_TestFunction_NoiseLevel);

	return fitness_value;
}

double Optimization_TestFunction_Schwefel2_26(double *solution)
{
	int i;
	double fitness_value;

	double sum_sinx;

	Optimization_TestFunction_RotationSolution(solution, 500.0, 1000.0);

	fitness_value = 0.0;
	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
		fitness_value -= Rotated_Solution[i] * sin(sqrt(fabs(Rotated_Solution[i])));

	if(Optimization_TestFunction_NoiseLevel >= 0.0 && Optimization_TestFunction_NoiseLevel <= 1000.0)
		fitness_value += randn(0, Optimization_TestFunction_NoiseLevel);


	sum_sinx = 0.0;
	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
		sum_sinx += (0.5 * (1.0 - cos(Rotated_Solution[i] * 1.0)));
	sum_sinx /= Optimization_TestFunction_Dimension;

	return fitness_value + 0.1 * sum_sinx * fitness_value;
}

double Optimization_TestFunction_Ackley(double *solution)
{
	int i;
	double sum_x, sum_cosx, fitness_value;

	double sum_sinx;

	Optimization_TestFunction_RotationSolution(solution, 32.0, 64.0);

	sum_x = 0.0;
	sum_cosx = 0.0;
	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
	{
		sum_x += pow(Rotated_Solution[i], 2);
		sum_cosx += cos(2.0 * PI * Rotated_Solution[i]);
	}

	fitness_value = -20.0 * exp(-0.2 * sqrt(sum_x / Optimization_TestFunction_Dimension)) -
					exp(sum_cosx / Optimization_TestFunction_Dimension) + 20.0 + exp(1);

//	if(Optimization_TestFunction_NoiseLevel >= 0.0 && Optimization_TestFunction_NoiseLevel <= 1000.0)
//		fitness_value += randn(0, Optimization_TestFunction_NoiseLevel);

	sum_sinx = 0.0;
	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
		sum_sinx += (0.5 * (1.0 - cos(Rotated_Solution[i] * 1.0)));
	sum_sinx /= Optimization_TestFunction_Dimension;

	return fitness_value + 0.1 * sum_sinx * fitness_value;
}

double Optimization_TestFunction_Foxholes(double *solution)
{
	int j;
	double x1, x2, fitness_value;

	if(Optimization_TestFunction_RotationFlag == 1)
	{
		Scaled_Solution[0] = 132.0 * solution[0] / Optimization_TestFunction_AxisRange - 98.0;
		Scaled_Solution[1] = 132.0 * solution[1] / Optimization_TestFunction_AxisRange - 98.0;

		x1 = Scaled_Solution[0] * Optimization_TestFunction_RotationMatrix[0][0] +
			 Scaled_Solution[1] * Optimization_TestFunction_RotationMatrix[1][0];
		x2 = Scaled_Solution[0] * Optimization_TestFunction_RotationMatrix[0][1] +
			 Scaled_Solution[1] * Optimization_TestFunction_RotationMatrix[1][1];
	}
	else
	{
		x1 = 132.0 * solution[0] / Optimization_TestFunction_AxisRange - 98.0;
		x2 = 132.0 * solution[1] / Optimization_TestFunction_AxisRange - 98.0;
	}

	fitness_value = 0.0;
	for(j=0; j< 25; ++j)
		fitness_value += 1.0 / ((j+1) +
								pow(x1 - Optimization_TestFunction_Foxholes_a1[j], 6) +
								pow(x2 - Optimization_TestFunction_Foxholes_a2[j], 6));

	fitness_value = 1.0 / (0.002 + fitness_value);

	if(Optimization_TestFunction_NoiseLevel >= 0.0 && Optimization_TestFunction_NoiseLevel <= 1000.0)
		fitness_value += randn(0, Optimization_TestFunction_NoiseLevel);

	return fitness_value;
}

double Optimization_TestFunction_SixHumpCamelBack(double *solution)
{
	double x1, x2, fitness_value;

	if(Optimization_TestFunction_RotationFlag == 1)
	{
		Scaled_Solution[0] = 10.0 * solution[0] / Optimization_TestFunction_AxisRange - 4.91017;
		Scaled_Solution[1] = 10.0 * solution[1] / Optimization_TestFunction_AxisRange - 5.7126;

		x1 = Scaled_Solution[0] * Optimization_TestFunction_RotationMatrix[0][0] +
			 Scaled_Solution[1] * Optimization_TestFunction_RotationMatrix[1][0];
		x2 = Scaled_Solution[0] * Optimization_TestFunction_RotationMatrix[0][1] +
			 Scaled_Solution[1] * Optimization_TestFunction_RotationMatrix[1][1];
	}
	else
	{
		x1 = 10.0 * solution[0] / Optimization_TestFunction_AxisRange - 4.91017;
		x2 = 10.0 * solution[1] / Optimization_TestFunction_AxisRange - 5.7126;
	}

	fitness_value = 4.0 * pow(x1, 2) -
					2.1 * pow(x1, 4) + 
					pow(x1, 6) / 3.0 +
					x1 * x2 -
					4.0 * pow(x2, 2) +
					4.0 * pow(x2, 4);

	if(Optimization_TestFunction_NoiseLevel >= 0.0 && Optimization_TestFunction_NoiseLevel <= 1000.0)
		fitness_value += randn(0, Optimization_TestFunction_NoiseLevel);

	return fitness_value;
}

double Optimization_TestFunction_Branin(double *solution)
{
	double x1, x2, fitness_value;

	if(Optimization_TestFunction_RotationFlag == 1)
	{
		Scaled_Solution[0] = 15.0 * solution[0] / Optimization_TestFunction_AxisRange - 8.142;
		Scaled_Solution[1] = 15.0 * solution[1] / Optimization_TestFunction_AxisRange - 12.275;

		x1 = Scaled_Solution[0] * Optimization_TestFunction_RotationMatrix[0][0] +
			 Scaled_Solution[1] * Optimization_TestFunction_RotationMatrix[1][0];
		x2 = Scaled_Solution[0] * Optimization_TestFunction_RotationMatrix[0][1] +
			 Scaled_Solution[1] * Optimization_TestFunction_RotationMatrix[1][1];
	}
	else
	{
		x1 = 15.0 * solution[0] / Optimization_TestFunction_AxisRange - 8.142;
		x2 = 15.0 * solution[1] / Optimization_TestFunction_AxisRange - 12.275;
	}

	fitness_value = pow(x2 -
						5 * pow(x1, 2) / (4.0 * pow(PI, 2)) +
						5.0 * x1 / PI -
						6.0, 2) +
					10.0 * (1.0 - 1.0 / (8.0 * PI)) * cos(x1) +
					10.0;

	if(Optimization_TestFunction_NoiseLevel >= 0.0 && Optimization_TestFunction_NoiseLevel <= 1000.0)
		fitness_value += randn(0, Optimization_TestFunction_NoiseLevel);

	return fitness_value;
}

double Optimization_TestFunction_GoldsteinPrice(double *solution)
{
	double x1, x2, g, h, fitness_value;

	if(Optimization_TestFunction_RotationFlag == 1)
	{
		Scaled_Solution[0] = 4.0 * solution[0] / Optimization_TestFunction_AxisRange - 2.0;
		Scaled_Solution[1] = 4.0 * solution[1] / Optimization_TestFunction_AxisRange - 2.0;

		x1 = Scaled_Solution[0] * Optimization_TestFunction_RotationMatrix[0][0] +
			 Scaled_Solution[1] * Optimization_TestFunction_RotationMatrix[1][0];
		x2 = Scaled_Solution[0] * Optimization_TestFunction_RotationMatrix[0][1] +
			 Scaled_Solution[1] * Optimization_TestFunction_RotationMatrix[1][1];
	}
	else
	{
		x1 = 4.0 * solution[0] / Optimization_TestFunction_AxisRange - 2.0;
		x2 = 4.0 * solution[1] / Optimization_TestFunction_AxisRange - 2.0;
	}

	g = 1.0 + pow(x1 + x2 +1.0, 2) * (19.0 - 14.0 * x1 + 3.0 * pow(x1, 2) - 14.0 * x2 + 6.0 * x1 * x2 + 3.0 * pow(x2, 2));
	h = 30.0 + pow(2.0 * x1 - 3.0 * x2, 2) * (18.0  - 32.0 * x1 + 12.0 * pow(x1, 2) + 48.0 * x2 - 36.0 * x1 * x2 + 27.0 * pow(x2, 2));
	
	fitness_value = g * h;

	if(Optimization_TestFunction_NoiseLevel >= 0.0 && Optimization_TestFunction_NoiseLevel <= 1000.0)
		fitness_value += randn(0, Optimization_TestFunction_NoiseLevel);

	return fitness_value;
}

double Optimization_TestFunction_HighConditionedElliptic(double *solution)
{
	int i;
	double fitness_value;

	Optimization_TestFunction_RotationSolution(solution, 100.0, 200.0);

	fitness_value = 0.0;
	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
		fitness_value += pow(10.0, 6 * (double) i / (Optimization_TestFunction_Dimension-1)) * pow(Rotated_Solution[i],2);

	if(Optimization_TestFunction_NoiseLevel >= 0.0 && Optimization_TestFunction_NoiseLevel <= 1.0)
		fitness_value += Optimization_TestFunction_NoiseLevel * fitness_value * rand() / RAND_MAX * 2.0 * (rand() / RAND_MAX - 0.5);

	return fitness_value;
}

double Optimization_TestFunction_Weierstrass(double *solution)
{
	int i, k;
	double a, b, fitness_value_1, fitness_value_2, fitness_value;

	Optimization_TestFunction_RotationSolution(solution, 0.5, 1.0);

	a = 0.5;
	b = 3.0;
	fitness_value_1 = 0.0;
	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
		for(k=0;k<=20; ++k)
			fitness_value_1 += (pow(a,k) * cos(2.0 * PI * pow(b,k) * (Rotated_Solution[i] + 0.5)));

//	fitness_value_2 = 0.0;
//	for(k=0;k<=20; ++k)
//		fitness_value_2 += (pow(a,k) * cos(2.0 * PI * pow(b,k) *  0.5));
	fitness_value_2 = -2.0;		//for k_max = 20;

	fitness_value = fitness_value_1 - Optimization_TestFunction_Dimension * fitness_value_2;

	if(Optimization_TestFunction_NoiseLevel >= 0.0 && Optimization_TestFunction_NoiseLevel <= 1.0)
		fitness_value += Optimization_TestFunction_NoiseLevel * fitness_value * rand() / RAND_MAX * 2.0 * (rand() / RAND_MAX - 0.5);

	return fitness_value;
}

double Optimization_TestFunction_Sinc(double *solution)
{
	int i;
	double fitness_value;

	Optimization_TestFunction_RotationSolution(solution, 1.0, 2.0);

	fitness_value = 0.0;
	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
		fitness_value += pow(Rotated_Solution[i], 2);

	fitness_value = sin(fitness_value * 7.0 * PI) / (fitness_value * 7.0 * PI);

	if(Optimization_TestFunction_NoiseLevel >= 0.0 && Optimization_TestFunction_NoiseLevel <= 1000.0)
		fitness_value += randn(0, Optimization_TestFunction_NoiseLevel);

	return fitness_value;
}

double Optimization_TestFunction_Levy(double *solution)
{
	int i;
	double fitness_value;

	Optimization_TestFunction_RotationSolution(solution, 10.0, 20.0);

	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
		Rotated_Solution[i] = 1.0 + (Rotated_Solution[i] - 1.0) / 4.0;

	fitness_value = pow(sin(PI * Rotated_Solution[0]), 2.0);

	for(i=0; i< Optimization_TestFunction_Dimension-1; ++i)
		fitness_value += (pow(Rotated_Solution[i] - 1.0, 2) *
						 (1.0 + 10.0 * pow(sin(PI * Rotated_Solution[i+1]), 2)));

	fitness_value += pow(Rotated_Solution[Optimization_TestFunction_Dimension-1] - 1.0, 2) *
					 (1.0 + 10.0 * pow(sin(2.0 * PI * Rotated_Solution[Optimization_TestFunction_Dimension-1]), 2));

	if(Optimization_TestFunction_NoiseLevel >= 0.0 && Optimization_TestFunction_NoiseLevel <= 1000.0)
		fitness_value += randn(0, Optimization_TestFunction_NoiseLevel);

	return fitness_value;
}

double Optimization_TestFunction_Pern(double *solution)
{
	int i, k;
	double cur_fitness_value, fitness_value;

	Optimization_TestFunction_RotationSolution(solution, Optimization_TestFunction_Dimension, 2.0*Optimization_TestFunction_Dimension);

	fitness_value = 0.0;
	for(k=0; k< Optimization_TestFunction_Dimension; ++k)
	{
		cur_fitness_value = 0.0;
		for(i=0; i< Optimization_TestFunction_Dimension; ++i)
			cur_fitness_value += (pow(i+1,(k+1)) + 0.5) * (pow((1.0 / (i+1) * Rotated_Solution[i]),k+1) - 1.0);
		fitness_value += pow(cur_fitness_value, 2);
	}

	if(Optimization_TestFunction_NoiseLevel >= 0.0 && Optimization_TestFunction_NoiseLevel <= 1000.0)
		fitness_value += randn(0, Optimization_TestFunction_NoiseLevel);

	return fitness_value;
}

double Optimization_TestFunction_Zakharov(double *solution)
{
	int i;
	double cur_fitness_value, fitness_value;

	Optimization_TestFunction_RotationSolution(solution, 5.0, 15.0);

	fitness_value = 0.0;
	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
			fitness_value += pow(Rotated_Solution[i], 2);

	cur_fitness_value = 0.0;
	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
		cur_fitness_value += 0.5 * (i + 1) * Rotated_Solution[i];
	fitness_value += pow(cur_fitness_value, 2) + pow(cur_fitness_value, 4);

	if(Optimization_TestFunction_NoiseLevel >= 0.0 && Optimization_TestFunction_NoiseLevel <= 1000.0)
		fitness_value += randn(0, Optimization_TestFunction_NoiseLevel);

	return fitness_value;
}

double Optimization_TestFunction_Alpine(double *solution)
{
	int i;
	double fitness_value;

	Optimization_TestFunction_RotationSolution(solution, 10.0, 20.0);

	fitness_value = 0.0;
	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
		fitness_value += fabs(Rotated_Solution[i] * sin(Rotated_Solution[i]) + 0.1 * Rotated_Solution[i]);

	return fitness_value;
}

double Optimization_TestFunction_Pathological(double *solution)
{
	int i;
	double fitness_value;

	Optimization_TestFunction_RotationSolution(solution, 100.0, 200.0);

	fitness_value = 0.0;
	for(i=0; i< Optimization_TestFunction_Dimension-1; ++i)
		fitness_value += 0.5 + (pow(sin(sqrt(100.0 * pow(Rotated_Solution[i],2) + pow(Rotated_Solution[i+1],2))), 2) - 0.5) /
							   (1.0 + 0.001 * pow(Rotated_Solution[i]-Rotated_Solution[i+1] ,4));

	return fitness_value;
}

double Optimization_TestFunction_InvertedCosineWave(double *solution)
{
	int i;
	double cur_fitness_value, fitness_value;

	Optimization_TestFunction_RotationSolution(solution, 5.0, 10.0);

	fitness_value = 0.0;
	for(i=0; i< Optimization_TestFunction_Dimension-1; ++i)
	{
		cur_fitness_value = pow(Rotated_Solution[i],2) + pow(Rotated_Solution[i+1],2) + 0.5 * Rotated_Solution[i] * Rotated_Solution[i+1];
		fitness_value -= exp(-cur_fitness_value / 8.0) *
						 cos(4.0 * sqrt(cur_fitness_value));
	}

	return fitness_value;
}

double Optimization_TestFunction_InvertedCosineMixture(double *solution)
{
	int i;
	double fitness_value;

	Optimization_TestFunction_RotationSolution(solution, 1.0, 2.0);

	fitness_value = 0.0;
	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
		fitness_value += 0.1 * cos(5.0 * PI * Rotated_Solution[i]);

	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
		fitness_value -= pow(Rotated_Solution[i], 2);

	fitness_value = 0.1 * Optimization_TestFunction_Dimension - fitness_value;

	return fitness_value;
}

double Optimization_TestFunction_EpistaticMichalewicz(double *solution)
{
	int i;
	double y, ph, m;
	double fitness_value;

	ph = PI / 6.0;
	m = 10;

	Optimization_TestFunction_RotationSolution(solution, 0.0, PI);

	fitness_value = 0.0;
	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
	{
		if(i == Optimization_TestFunction_Dimension - 1)
			y = Rotated_Solution[i];
		else if((i+1) % 2 == 1)
			y = Rotated_Solution[i] * cos(ph) - Rotated_Solution[i+1] * sin(ph);
		else
			y = Rotated_Solution[i] * sin(ph) + Rotated_Solution[i+1] * cos(ph);

		fitness_value -= sin(y) * pow(sin((i+1) * pow(y,2) / PI), 2.0 * m);
	}

	return fitness_value;
}

double Optimization_TestFunction_LevyMontalvo(double *solution)
{
	int i;
	double fitness_value;

	Optimization_TestFunction_RotationSolution(solution, 5.0, 10.0);

	fitness_value = pow(sin(3.0 * PI * Rotated_Solution[0]), 2);
	for(i=0; i< Optimization_TestFunction_Dimension-1; ++i)
		fitness_value += pow(Rotated_Solution[i] - 1.0, 2) * (1.0 + pow(3.0 * PI * Rotated_Solution[i+1], 2));

	fitness_value += pow(Rotated_Solution[Optimization_TestFunction_Dimension-1] - 1.0, 2) * (1.0 + pow(2.0 * PI * Rotated_Solution[Optimization_TestFunction_Dimension-1], 2));

	return fitness_value;
}

double Optimization_TestFunction_Neumaier3(double *solution)
{
	int i;
	double fitness_value;

	Optimization_TestFunction_RotationSolution(solution, pow(Optimization_TestFunction_Dimension,2), 2.0*pow(Optimization_TestFunction_Dimension,2));

//	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
//		Rotated_Solution[i] = (i+1)*(Optimization_TestFunction_Dimension+1-(i+1));

	fitness_value = 0.0;
	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
		fitness_value += pow(Rotated_Solution[i] - 1.0, 2);

	for(i=1; i< Optimization_TestFunction_Dimension; ++i)
		fitness_value -= Rotated_Solution[i] * Rotated_Solution[i-1];

//	printf("%lf ", fitness_value);
//	Pause();
	return fitness_value;
}

double Optimization_TestFunction_OddSquare(double *solution)
{
	int i;
	double d, D, fitness_value;

	Optimization_TestFunction_RotationSolution(solution, 15.0, 30.0);

//	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
//		Rotated_Solution[i] = Optimization_TestFunction_OddSquare_b[i % 10];

	d = 0.0;
	D = 0.0;
	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
	{
		d += pow(Rotated_Solution[i] - Optimization_TestFunction_OddSquare_b[i % 10], 2);
		if(D < fabs(Rotated_Solution[i] - Optimization_TestFunction_OddSquare_b[i % 10]))
			D = fabs(Rotated_Solution[i] - Optimization_TestFunction_OddSquare_b[i % 10]);
	}
	d = sqrt(d);
	D *= sqrt(Optimization_TestFunction_Dimension);

	fitness_value = -(1.0 + 0.2*d / (D + 0.1)) * cos(D * PI) * exp(-D / (2.0 * PI));

//	printf("%lf ", fitness_value);
//	Pause();
	return fitness_value;
}

double Optimization_TestFunction_Paviani(double *solution)
{
	int i;
	double cur_fitness_value, fitness_value;

	Optimization_TestFunction_RotationSolution(solution, -2.0, 8.0);

	fitness_value = 0.0;
	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
		fitness_value += pow(log(Rotated_Solution[i] - 2.0), 2) + pow(log(10.0 - Rotated_Solution[i]), 2);

	cur_fitness_value = 1.0;
	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
		cur_fitness_value *= Rotated_Solution[i];

	fitness_value = fitness_value - pow(cur_fitness_value, 0.2);

	return fitness_value;
}

double Optimization_TestFunction_Periodic(double *solution)
{
	int i;
	double cur_fitness_value, fitness_value;

	Optimization_TestFunction_RotationSolution(solution, 10.0, 20.0);

	fitness_value = 1.0;
	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
		fitness_value += pow(sin(Rotated_Solution[i]), 2);

	cur_fitness_value = 0.0;
	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
		cur_fitness_value -= pow(Rotated_Solution[i], 2);

	fitness_value -= 0.1 * exp(cur_fitness_value);

	return fitness_value;
}

double Optimization_TestFunction_Salomon(double *solution)
{
	int i;
	double fitness_value;

	Optimization_TestFunction_RotationSolution(solution, 100.0, 200.0);

	fitness_value = 0.0;
	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
		fitness_value += pow(Rotated_Solution[i], 2);
	fitness_value = sqrt(fitness_value);

	fitness_value = 1.0 - cos(2.0 * PI * fitness_value) + 0.1 * fitness_value;

	return fitness_value;
}

double Optimization_TestFunction_Shubert(double *solution)
{
	int i, j;
	double cur_fitness_value, fitness_value;

	Optimization_TestFunction_RotationSolution(solution, 10.0, 20.0);

	fitness_value = 1.0;
	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
	{
		cur_fitness_value = 0.0;
		for(j=0; j< 5; ++j)
			cur_fitness_value += (j + 1.0) * cos((j + 2.0) * Rotated_Solution[i] + (j + 1.0));
		fitness_value *= cur_fitness_value;
	}

	return fitness_value;
}

double Optimization_TestFunction_Sinusoidal(double *solution)
{
	int i;
	double A, B, z;
	double f1, f2, fitness_value;

	A = 2.5;
	B = 5.0;
	z = 30.0;

	Optimization_TestFunction_RotationSolution(solution, 0.0, 180.0);

	f1 = f2 = 1.0;
	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
	{
		f1 *= sin((Rotated_Solution[i] - z) * PI / 180.0);
		f2 *= sin(B * (Rotated_Solution[i] - z) * PI / 180.0);
	}

	fitness_value = -(A*f1 + f2);

	return fitness_value;
}

double Optimization_TestFunction_Michalewicz(double *solution)
{
	int i;
	double m, fitness_value;

	Optimization_TestFunction_RotationSolution(solution, 0.0, PI);

	m = 10.0;

	fitness_value = 0.0;
	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
		fitness_value -= sin(Rotated_Solution[i]) * pow(sin((i+1.0) * pow(Rotated_Solution[i], 2) / PI), 2.0 * m);

	return fitness_value;
}

double Optimization_TestFunction_Whitely(double *solution)
{
	int i, j;
	double one_x2, y_ij, fitness_value;

	Optimization_TestFunction_RotationSolution(solution, 100.0, 200.0);

	fitness_value = 0.0;
	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
	{
		one_x2 = pow(1.0 - Rotated_Solution[i], 2);
		for(j=0; j< Optimization_TestFunction_Dimension; ++j)
		{
			y_ij = 100.0 * pow(Rotated_Solution[j] - pow(Rotated_Solution[i], 2), 2) + one_x2;
			fitness_value += y_ij / 4000 - cos(y_ij) + 1.0;
		}
	}

	return fitness_value;
}

double Optimization_TestFunction_Easom(double *solution)
{
	int i;
	double sum_x_pi, fitness_value;

	Optimization_TestFunction_RotationSolution(solution, 10.0, 20.0);

	fitness_value = 1.0;
	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
		fitness_value *= cos(Rotated_Solution[i]);

	sum_x_pi = 0.0;
	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
		sum_x_pi -= pow(Rotated_Solution[i] - PI, 2);

	fitness_value = -fitness_value * exp(sum_x_pi);

	return fitness_value;
}

//================================================================================================================
//==================================== Hybrid Composition Objective Function =====================================
//================================================================================================================

double Optimization_TestFunction_HybridComposition01(double *solution)
{
	int i, j, k;
	double a, b, C;
	double x, sum_x, sum_cosx, product_cosx; 
	double *cur_solution;
	double w[5], max_w, sum_w, w_weight;
	double cur_f, cur_f1, cur_f2, f[10], f_max;
	double fitness_value;

	cur_solution = (double *)malloc(sizeof(double) * Optimization_TestFunction_Dimension);

	C = 2000.0;

	for(i=0; i< 5; ++i)
	{
		w[i] = 0.0;
		for(j=0; j< Optimization_TestFunction_Dimension; ++j)
			w[i] += pow((10.0*solution[j] / Optimization_TestFunction_AxisRange-5.0) - i*0.5, 2);
		w[i] = exp(-1.0 * w[i] / (2 * Optimization_TestFunction_Dimension * pow(1, 2)));

		if(i==0 || max_w < w[i])
			max_w = w[i];
	}

	sum_w = 0.0;
	w_weight = 1.0 - pow(max_w, 10);
	for(i=0; i< 5; ++i)
	{
		if(w[i] < max_w)
			w[i] = w[i] * w_weight;

		sum_w += w[i];
	}

	for(i=0; i< 5; ++i)
		w[i] /= sum_w;

	//------------ f1 - f2
	cur_f = 0.0;
	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
	{
		x = 10.0 * solution[i] / Optimization_TestFunction_AxisRange - 5.0 + 0.0;
		cur_f += (pow(x, 2) -
			      10.0 * cos(2.0 * PI * x) +
			 	 10.0);
	}

	f_max = 0.0;
	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
	{
		x = 5.0 + 0.0;
		f_max += (pow(x, 2) -
			      10.0 * cos(2.0 * PI * 5.0) +
			 	 10.0);
	}

	f[0] = w[0] * (C*cur_f/f_max + 0.0);
	f[1] = w[0] * (C*cur_f/f_max + 100.0);

	//------------ f3 - f4
	a = 0.5;
	b = 3.0;
	cur_f1 = 0.0;
	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
	{
		x = 10.0 * (10.0 * solution[i] / Optimization_TestFunction_AxisRange - 5.0 + 0.5);

		for(k=0;k<=20; ++k)
			cur_f1 += (pow(a,k) * cos(2.0 * PI * pow(b,k) * (x + 0.5)));
	}

	cur_f2 = 0.0;
	for(k=0;k<=20; ++k)
		cur_f2 += (pow(a,k) * cos(2.0 * PI * pow(b,k) *  0.5));

	cur_f = cur_f1 - Optimization_TestFunction_Dimension * cur_f2;

	cur_f1 = 0.0;
	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
	{
		x = 10.0 * (5.0 + 0.5);

		for(k=0;k<=20; ++k)
			cur_f1 += (pow(a,k) * cos(2.0 * PI * pow(b,k) * (x + 0.5)));
	}

	cur_f2 = 0.0;
	for(k=0;k<=20; ++k)
		cur_f2 += (pow(a,k) * cos(2.0 * PI * pow(b,k) *  0.5));

	f_max = cur_f1 - Optimization_TestFunction_Dimension * cur_f2;

	f[2] = w[1] * (C*cur_f/f_max + 200.0);
	f[3] = w[1] * (C*cur_f/f_max + 300.0);

	//------------ f5 - f6
	sum_x = 0.0;
	product_cosx = 1.0;
	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
	{
		x = 12.0 * (10.0 * solution[i] / Optimization_TestFunction_AxisRange - 5.0 + 1.0);
		sum_x += pow(x, 2);
		product_cosx *= cos(x / sqrt(i+1));
	}

	cur_f = sum_x/4000.0 - product_cosx + 1.0;

	sum_x = 0.0;
	product_cosx = 1.0;
	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
	{
		x = 12.0 * (5.0 + 1.0);
		sum_x += pow(x, 2);
		product_cosx *= cos(x / sqrt(i+1));
	}

	f_max = sum_x/4000.0 - product_cosx + 1.0;

	f[4] = w[2] * (C*cur_f/f_max + 400.0);
	f[5] = w[2] * (C*cur_f/f_max + 500.0);

	//------------ f7 - f8
	sum_x = 0.0;
	sum_cosx = 0.0;
	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
	{
		x = 6.4 * (10.0 * solution[i] / Optimization_TestFunction_AxisRange - 5.0 + 1.5);
		sum_x += pow(x, 2);
		sum_cosx += cos(2.0 * PI * x);
	}

	cur_f = -20.0 * exp(-0.2 * sqrt(sum_x / 30.0)) -
			exp(sum_cosx / 30.0) + 20.0 + exp(1);

	sum_x = 0.0;
	sum_cosx = 0.0;
	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
	{
		x = 6.4 * (5.0 + 1.5);

		sum_x += pow(x, 2);
		sum_cosx += cos(2.0 * PI * x);
	}

	f_max = -20.0 * exp(-0.2 * sqrt(sum_x / 30.0)) -
			exp(sum_cosx / 30.0) + 20.0 + exp(1);

	f[6] = w[3] * (C*cur_f/f_max + 600.0);
	f[7] = w[3] * (C*cur_f/f_max + 700.0);

	//------------ f9 - f10
	cur_f = 0.0;
	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
	{
		x = 20.0 * (10.0 * solution[i] / Optimization_TestFunction_AxisRange - 5.0 + 2.0);
		cur_f += pow(x, 2);
	}

	f_max = 0.0;
	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
	{
		x = 20.0 * (5.0 + 2.0);
		f_max += pow(x, 2);
	}

	f[8] = w[4] * (C*cur_f/f_max + 800.0);
	f[9] = w[4] * (C*cur_f/f_max + 900.0);

//	for(i=0; i<10; ++i)
//		printf("%1.9lf ", f[i]);
//	Pause();

	free(cur_solution);

	fitness_value = 0.0;
	for(i=0; i < 10; ++i)
		fitness_value += f[i];

	if(Optimization_TestFunction_NoiseLevel >= 0.0 && Optimization_TestFunction_NoiseLevel <= 1000.0)
		fitness_value += randn(0, Optimization_TestFunction_NoiseLevel);

	return fitness_value;
}

//================================================================================================================
//==================================== Pseudo Boolean Objective Function =========================================
//================================================================================================================

double Optimization_TestFunction_MaxOne(double *solution)
{
	int i;
	double fitness_value;

	fitness_value = 0.0;
	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
	{
		if(solution[i] > 0.5)
			++fitness_value;
	}

	return fitness_value;
}

double Optimization_TestFunction_RealRoyalRoad(double *solution)
{
	int i;
	double no_one, fitness_value;

	if(Optimization_TestFunction_RealRoyalRoadTau < 1 || Optimization_TestFunction_RealRoyalRoadTau > Optimization_TestFunction_Dimension)
		Optimization_TestFunction_RealRoyalRoadTau = (int) (Optimization_TestFunction_Dimension / 2);

	no_one = 0.0;
	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
		if(solution[i] >= 0.5)
			++no_one;

	if(no_one == Optimization_TestFunction_Dimension || no_one <= Optimization_TestFunction_Dimension - Optimization_TestFunction_RealRoyalRoadTau)
		fitness_value = Optimization_TestFunction_RealRoyalRoadTau + no_one;
	else
		fitness_value = Optimization_TestFunction_Dimension - no_one;

//	printf("%d %d %lf %lf\n", Optimization_TestFunction_Dimension, Optimization_TestFunction_RealRoyalRoadTau, no_one, fitness_value);
//	scanf("%d", &i);

	return Optimization_TestFunction_Dimension - fitness_value;
}

double Optimization_TestFunction_InverseRealRoyalRoad2001(double *solution)
{
	int i;
	int no_one, cur_block_length, max_block_length;

	if(Optimization_TestFunction_RealRoyalRoadTau < 1 || Optimization_TestFunction_RealRoyalRoadTau > Optimization_TestFunction_Dimension)
		Optimization_TestFunction_RealRoyalRoadTau = (int) (Optimization_TestFunction_Dimension / 2);

	no_one = 0;
	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
		no_one += (int) solution[i];

	if(no_one == Optimization_TestFunction_Dimension)
		return -2.0 * pow(Optimization_TestFunction_Dimension, 2);
	else
	{
		max_block_length = 0;
		cur_block_length = 0;

		for(i=0; i< Optimization_TestFunction_Dimension; ++i)
		{
			if(solution[i] < 0.5)
			{
				if(max_block_length < cur_block_length)
					max_block_length = cur_block_length;
				cur_block_length = 0;
			}
			else
				++cur_block_length;
		}

		if(max_block_length < cur_block_length)
			max_block_length = cur_block_length;

		if(no_one <= Optimization_TestFunction_Dimension - Optimization_TestFunction_RealRoyalRoadTau)
			return -2 * no_one - max_block_length;

		return 0;
	}

	return 0;
}

double Optimization_TestFunction_Linear(double *solution)
{
	int i;
	double fitness_value;

//	fitness_value = Optimization_TestFunction_LinearWeigth[0];
//	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
//		fitness_value += Optimization_TestFunction_LinearWeigth[i+1] * solution[i];

	fitness_value = 0;
	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
		fitness_value += (i+1) * solution[i];

	return -1*fitness_value/(0.5*Optimization_TestFunction_Dimension*(Optimization_TestFunction_Dimension+1));
}

double Optimization_TestFunction_AlmostPositive(double *solution)
{
	int i;
	double sum_gene, prd_gene, fitness_value;

	sum_gene = 0;
	prd_gene = 1;
	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
	{
		sum_gene += solution[i];
		prd_gene *= solution[i];
	}

	fitness_value = Optimization_TestFunction_Dimension - sum_gene +
					(Optimization_TestFunction_Dimension + 1) * prd_gene;

	return fitness_value;
}

double Optimization_TestFunction_InversePseudoModular(double *solution)
{
	int i;
	double fitness_value;

	fitness_value = 0;
	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
		if(solution[i] == 1.0)
			++fitness_value;
		else
			break;

	return Optimization_TestFunction_Dimension - fitness_value;
}

double Optimization_TestFunction_RealValue(double *solution)
{
	int i;
	double fitness_value;

	fitness_value = 0;
	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
		if(solution[i] == 1.0)
			fitness_value += pow(2, i);

	return fabs(pow(2, Optimization_TestFunction_Dimension) - fitness_value);
}

double Optimization_TestFunction_MultiRealValue(double *solution)
{
	int i;
	int min_gene_index, max_gene_index;
	double cur_fitness_value, fitness_value;

	if(Optimization_TestFunction_MultiRealValue_BlockSize < 1 || Optimization_TestFunction_MultiRealValue_BlockSize > Optimization_TestFunction_Dimension)
		Optimization_TestFunction_MultiRealValue_BlockSize = (int) (Optimization_TestFunction_Dimension / 2);

	min_gene_index = 0;
	max_gene_index = Optimization_TestFunction_MultiRealValue_BlockSize;

	fitness_value = 0;
	while (1) {
		cur_fitness_value = 0;
		for(i=min_gene_index; i< max_gene_index; ++i)
			if(solution[i] == 1.0)
				cur_fitness_value += pow(2, i-min_gene_index);
		fitness_value += fabs(cur_fitness_value - pow(2, max_gene_index - min_gene_index));
		
		if(max_gene_index == Optimization_TestFunction_Dimension)
			break;

		min_gene_index += Optimization_TestFunction_MultiRealValue_BlockSize;
		max_gene_index += Optimization_TestFunction_MultiRealValue_BlockSize;
		if(max_gene_index >= Optimization_TestFunction_Dimension)
			max_gene_index = Optimization_TestFunction_Dimension;
	}

	return fitness_value;
}

double Optimization_TestFunction_ShortPathConstant(double *solution)
{
	int i, mode;
	double sum_One, fitness_value;

	sum_One = 0.0;
	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
	{
		if(solution[i] > 0.5)
			++sum_One;
	}

	if(sum_One == Optimization_TestFunction_Dimension)
		fitness_value = 2.0 * Optimization_TestFunction_Dimension;
	else
	{
		mode = 0;
		for(i=0; i< Optimization_TestFunction_Dimension - 1; ++i)
			if(solution[i] < 0.5 && solution[i+1] > 0.5)
			{
				mode = 1;
				break;
			}

		if(mode == 0)
			fitness_value = Optimization_TestFunction_Dimension + 1;
		else
			fitness_value = Optimization_TestFunction_Dimension - sum_One;
	}

	return Optimization_TestFunction_Dimension - fitness_value;
}

double Optimization_TestFunction_XORTree(double *solution)
{
	int i, cur_no_dimension, nxt_no_dimension;
	double fitness_value;

	fitness_value = 0;

	cur_no_dimension = Optimization_TestFunction_Dimension + 1;
	Optimization_TestFunction_XORTReeNode[0][0] = 1;
	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
		Optimization_TestFunction_XORTReeNode[0][i+1] = solution[i];

//	for(i=0; i< cur_no_dimension; ++i)
//		printf("%d ", (int) Optimization_TestFunction_XORTReeNode[0][i]);
//	printf("\n");

	do {
		for(i=0; i< (int) (cur_no_dimension / 2); ++i)
		{
		  if(Optimization_TestFunction_XORTReeNode[0][2*i] == Optimization_TestFunction_XORTReeNode[0][2*i+1])
		  {
			  Optimization_TestFunction_XORTReeNode[1][i] = 1;
			  --fitness_value;
		  }
		  else
			  Optimization_TestFunction_XORTReeNode[1][i] = 0;
		}
		nxt_no_dimension = (int) (cur_no_dimension / 2);

		if(cur_no_dimension % 2 == 1)
		{
			Optimization_TestFunction_XORTReeNode[1][nxt_no_dimension] = Optimization_TestFunction_XORTReeNode[0][cur_no_dimension-1];
			++nxt_no_dimension;
		}

		cur_no_dimension = nxt_no_dimension;
		for(i=0; i< cur_no_dimension; ++i)
			Optimization_TestFunction_XORTReeNode[0][i] = Optimization_TestFunction_XORTReeNode[1][i];

//		for(i=0; i< cur_no_dimension; ++i)
//			printf("%d ", (int) Optimization_TestFunction_XORTReeNode[0][i]);
//		printf("\n");

	} while(cur_no_dimension > 1);

//	printf("%d\n", (int) fitness_value);
//	scanf("%d", &i);

	return fitness_value * -1.0;
}


double Optimization_TestFunction_TruthTable(double *solution)
{
	int i;
	double dec_index;

	dec_index = 0;
	for(i=0; i < Optimization_TestFunction_Dimension; ++i)
		if(solution[i] == 1)
			dec_index += pow(2, i);

	return Optimization_TestFunction_TruthTableValue[(int) dec_index];
}

void Optimization_TestFunction_TruthTableLoad(char *truthtable_filename)
{
	int i;
	FILE *truthtable_ptr;

	truthtable_ptr = fopen(truthtable_filename, "r");
	
	fscanf(truthtable_ptr, "%d", &Optimization_TestFunction_Dimension);

	Optimization_TestFunction_TruthTableValue = (double *)malloc(sizeof(double) * (int) pow(2, Optimization_TestFunction_Dimension));

	for(i = 0; i < (int) pow(2, Optimization_TestFunction_Dimension); ++i)
		fscanf(truthtable_ptr, "%lf", &Optimization_TestFunction_TruthTableValue[i]);

	fclose(truthtable_ptr);

	return;
}

void Optimization_TestFunction_TruthTableConvert(double (*fn_ptr)(double*), int no_dimension)
{
	int dec_index;
	double *solution;
	
	solution = (double *)malloc(sizeof(double) * no_dimension);

	Optimization_TestFunction_Dimension = no_dimension;

	Optimization_TestFunction_TruthTableValue = (double *)malloc(sizeof(double) * (int) pow(2, no_dimension));

	for(dec_index = 0; dec_index < (int) pow(2, no_dimension); ++dec_index)
	{
		Dec2Bin(dec_index, &solution, no_dimension);
		Optimization_TestFunction_TruthTableValue[dec_index] = (*fn_ptr)(solution);
	}

	free(solution);

	return;
}

void Optimization_TestFunction_TruthTableRandomSwapping()
{
	int i, j;
	int **solution_indextree, *no_solution, cur_hamming;
	int solution_index1, solution_index2, solution_indextmp;
	double *solution, tmp_fitness;

	solution = (double *)malloc(sizeof(double) * Optimization_TestFunction_Dimension);
	no_solution = (int *)malloc(sizeof(int) * (Optimization_TestFunction_Dimension + 1));
	solution_indextree = (int **)malloc(sizeof(int*) * (Optimization_TestFunction_Dimension + 1));
	for(i=0; i< Optimization_TestFunction_Dimension + 1; ++i)
	{
		no_solution[i] = 0;
		solution_indextree[i] = (int *)malloc(sizeof(int) * (int) nChoosek(Optimization_TestFunction_Dimension, i));
	}

	for(i=0; i<(int) pow(2, Optimization_TestFunction_Dimension); ++i)
	{
		Dec2Bin(i, &solution, Optimization_TestFunction_Dimension);
		cur_hamming = 0;
		for(j=0; j< Optimization_TestFunction_Dimension; ++j)
			cur_hamming += (int) solution[j];

		solution_indextree[cur_hamming][no_solution[cur_hamming]] = i;
		++no_solution[cur_hamming];
	}

	for(i=1; i< Optimization_TestFunction_Dimension-1; ++i)
	{
		for(j=0; j<no_solution[i] / 2; ++j)
		{
			solution_index1 = rand() % no_solution[i];
			solution_index2 = rand() % no_solution[i];
			while(solution_index1 == solution_index2)
				solution_index2 = rand() % no_solution[i];

			solution_indextmp = solution_indextree[i][solution_index1];
			solution_indextree[i][solution_index1] = solution_indextree[i][solution_index2];
			solution_indextree[i][solution_index2] = solution_indextmp;

			tmp_fitness = Optimization_TestFunction_TruthTableValue[solution_indextree[i][solution_index1]];
			Optimization_TestFunction_TruthTableValue[solution_indextree[i][solution_index1]] = Optimization_TestFunction_TruthTableValue[solution_indextree[i][solution_index2]];
			Optimization_TestFunction_TruthTableValue[solution_indextree[i][solution_index2]] = tmp_fitness;
		}
	}

	for(i=0; i< Optimization_TestFunction_Dimension + 1; ++i)
		free(solution_indextree[i]);
	free(solution_indextree);
	free(no_solution);
	free(solution);

	return;
}

void Optimization_TestFunction_TruthTableDestruction()
{
	free(Optimization_TestFunction_TruthTableValue);

	return;
}

//================================================================================================================
//==================================== Fitness Landscape Rotation ================================================
//================================================================================================================

void Optimization_TestFunction_RotationMatrix_Construction()
{
	int i;

	Optimization_TestFunction_RotationFlag = 0;

	Scaled_Solution = (double *)malloc(sizeof(double) * Optimization_TestFunction_Dimension);
	Rotated_Solution = (double *)malloc(sizeof(double) * Optimization_TestFunction_Dimension);

	Optimization_TestFunction_RotationMatrix = (double **)malloc(sizeof(double*) * Optimization_TestFunction_Dimension);
	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
		Optimization_TestFunction_RotationMatrix[i] = (double *)malloc(sizeof(double) * Optimization_TestFunction_Dimension);

	OrthogonalMatrix_Generator(Optimization_TestFunction_Dimension, &Optimization_TestFunction_RotationMatrix);

	return;
}

void Optimization_TestFunction_RotationMatrix_Load()
{
	char rotate_filename[256];
	int i, j;
	FILE *rotate_fileptr;

	Optimization_TestFunction_RotationFlag = 1;

	Scaled_Solution = (double *)malloc(sizeof(double) * Optimization_TestFunction_Dimension);
	Rotated_Solution = (double *)malloc(sizeof(double) * Optimization_TestFunction_Dimension);

	Optimization_TestFunction_RotationMatrix = (double **)malloc(sizeof(double*) * Optimization_TestFunction_Dimension);
	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
		Optimization_TestFunction_RotationMatrix[i] = (double *)malloc(sizeof(double) * Optimization_TestFunction_Dimension);

	if(Optimization_TestFunction_Dimension < 10)
		sprintf(rotate_filename, "M_00%d", Optimization_TestFunction_Dimension);
	else
		sprintf(rotate_filename, "M_0%d", Optimization_TestFunction_Dimension);

	rotate_fileptr = fopen(rotate_filename, "r");
	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
		for(j=0; j< Optimization_TestFunction_Dimension; ++j)
			fscanf(rotate_fileptr, "%lf", &Optimization_TestFunction_RotationMatrix[i][j]);

	fclose(rotate_fileptr);

	return;
}

void Optimization_TestFunction_RotationMatrix_Destruction()
{
	int i;

	Optimization_TestFunction_RotationFlag = 0;

	free(Scaled_Solution);
	free(Rotated_Solution);
	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
		free(Optimization_TestFunction_RotationMatrix[i]);
	free(Optimization_TestFunction_RotationMatrix);

	return;
}

void Optimization_TestFunction_RotationSolution(double *cur_Solution, double min_x, double dx)
{
	int i, j;

	for(i=0; i< Optimization_TestFunction_Dimension; ++i)
		Scaled_Solution[i] = dx * cur_Solution[i] / Optimization_TestFunction_AxisRange - min_x;

	if(Optimization_TestFunction_RotationFlag == 1)
	{
		for(i=0; i< Optimization_TestFunction_Dimension; ++i)
		{
			Rotated_Solution[i] = 0.0;
			for(j=0; j< Optimization_TestFunction_Dimension; ++j)
				Rotated_Solution[i] += Scaled_Solution[j] * Optimization_TestFunction_RotationMatrix[j][i];
		}
	}
	else
		for(i=0; i< Optimization_TestFunction_Dimension; ++i)
			Rotated_Solution[i] = Scaled_Solution[i];

	return;
}

