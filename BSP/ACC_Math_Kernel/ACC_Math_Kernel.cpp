/******************************************************************************************/
//	Program Name	:	Non-Revisiting Genetic Algorithm (NrGA)
//	File Name		:	ACC_Math_Kernel.c
//	Auther			:	Dr. Chow Chi Kin, Dr. Yuen Shiu Yin
//	Edit by			:	Leung Shing Wa
//	University		:	City University of Hong Kong
//	Department		:	Electronic Engineering
//	Last Update		:	10 Sep 2009
//	Reference		:	A Genetic Algorithm that Adaptively Mutates and Never Revisits, 
//						IEEE Transactions on Evolutionary Computation, 
//						Vol 13(2) (April 2009) 454-472.
//
//	Discription		:	Functions that require in calculations are defined in this file.
#include "../../BSP/ACC_Math_Kernel/ACC_Math_Kernel.h"

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<malloc.h>
#include<time.h>

double Exponential_Width_neg, Exponential_Width_pos, *Exponential_Array;
double Normalize_Square_Width, *Normalize_Square_Array;

void Exponential_Construction()
{
	long l;

	Exponential_Array = (double *)malloc(sizeof(double) * EXPONENTIAL_NO_DIVID);
	Exponential_Width_neg = -1.0 * pow(3.0, 2) / (EXPONENTIAL_NO_DIVID - 1);

	for(l=0; l< EXPONENTIAL_NO_DIVID; ++l)
		Exponential_Array[l] = exp((double) l * Exponential_Width_neg);

	Exponential_Width_pos = -1.0 * Exponential_Width_neg;

	return;
}

void Exponential_Destruction()
{
	free(Exponential_Array);

	return;
}

double Exponential(double x)
{
	double exponential_index;

	if(x <=0)
		exponential_index = x / Exponential_Width_neg;
	else
		exponential_index = x / Exponential_Width_pos;
	
	if(exponential_index >= EXPONENTIAL_NO_DIVID)
		return 0.0;
	return Exponential_Array[(long) (exponential_index)];
}

void Normalize_Square_Construction()
{
	long l;

	Normalize_Square_Array = (double *)malloc(sizeof(double) * NORMALIZE_SQUARE_NO_DIVID);
	Normalize_Square_Width = 1.0 / (NORMALIZE_SQUARE_NO_DIVID - 1);

	for(l=0; l< NORMALIZE_SQUARE_NO_DIVID; ++l)
		Normalize_Square_Array[l] = pow((double) l * Normalize_Square_Width, 2);

	return;
}

void Normalize_Square_Destruction()
{
	free(Normalize_Square_Array);

	return;
}

double Normalize_Square(double x)
{
	long normalize_square_index;

	normalize_square_index = (long) (x / Normalize_Square_Width);

	return Normalize_Square_Array[abs(normalize_square_index)];
}

double randu()
{
	int i;
	double randu_value;

	randu_value = 0.0;
	for(i=0; i< 10; ++i)
		randu_value += rand();

	return randu_value / (RAND_MAX * 10.0);
}

double randn(double mean, double stdev)
{
	int i;
	double randn_value = 0;

	for(i=0; i< 12; ++i)
		randn_value += ((double) rand() / RAND_MAX - 0.5);
	
	return mean + stdev * randn_value;
}

void Sorting_Indexed(double **value, int **value_index, int no_item)
{
	int i, j, k;
	int sedgewick[] ={929, 505, 209, 41, 19, 5, 1}, increment, tmp_index;
	double tmp_value;

	for(k = 0; k < 7; ++k)
	{
		increment = sedgewick[k];

		for(i = increment; i < no_item; ++i)
		{
			tmp_value = *(*value + i);
			tmp_index = *(*value_index + i);
			for(j = i; j >= increment; j -= increment)
			{
				if (tmp_value < *(*value + j - increment))
				{
					*(*value + j) = *(*value + j - increment);
					*(*value_index + j)  = *(*value_index + j - increment);
				}
				else
					break;
			}
			*(*value + j) = tmp_value;
			*(*value_index + j) = tmp_index;
		}
	}

	return;
}

void Sorting_UnIndexed(double **value, int no_item)
{
	int i, j, k;
	int sedgewick[] ={929, 505, 209, 41, 19, 5, 1}, increment;
	double tmp_value;

	for(k = 0; k < 7; ++k)
	{
		increment = sedgewick[k];

		for(i = increment; i < no_item; ++i)
		{
			tmp_value = *(*value + i);
			for(j = i; j >= increment; j -= increment)
			{
				if (tmp_value < *(*value + j - increment))
					*(*value + j) = *(*value + j - increment);
				else
					break;
			}
			*(*value + j) = tmp_value;
		}
	}

	return;
}


double IMGau_Gau(double *Mu_1, double *Mu_2, double Var_1, double Var_2, int no_dimension, double **range)
{
	int k;
	double IGau_Gau_val = 1.0;

	for(k = 0; k < no_dimension; ++k)
		IGau_Gau_val *= IGau_Gau(Mu_1[k], Mu_2[k], Var_1, Var_2, range[k][0], range[k][1]);

	return IGau_Gau_val;
}

double IMGau(double *Mu, double Var, int no_dimension, double **range)
{
	int k;
	double IGau_val = 1.0;

	for(k=0; k< no_dimension; ++k)
		IGau_val *= IGau(Mu[k], Var, range[k][0], range[k][1]);

	return IGau_val;
}

double IGau_Gau(double Mu_1, double Mu_2, double Var_1, double Var_2, double Min_Range, double Max_Range)
{
	double bt, ap, M;
	double IGau_Gau_val;


	bt = (pow(Var_1, 2) * Mu_2 + pow(Var_2, 2) * Mu_1) / (pow(Var_1, 2) + pow(Var_2, 2));
//	bt = (var_1^2 * mu_2 + var_2^2 * mu_1) / (var_1^2 + var_2^2);
	ap = Var_1 * Var_2 / sqrt(pow(Var_1, 2) + pow(Var_2, 2));
//	ap = var_1 * var_2 / sqrt(var_1^2 + var_2^2);
	M = exp(-1.0 * ((pow(Var_1 * Mu_2, 2) + pow(Var_2 * Mu_1, 2)) / (pow(Var_1, 2) + pow(Var_2, 2)) - pow((pow(Var_1, 2) * Mu_2 + pow(Var_2, 2) * Mu_1) / (pow(Var_1, 2) + pow(Var_2, 2)), 2)) / (2.0 * pow(ap, 2)));
//	M = exp(-1 * ((var_1^2 * mu_2^2 + var_2^2 * mu_1^2) / (var_1^2 + var_2^2) - ((var_1^2 * mu_2 + var_2^2 * mu_1) / (var_1^2 + var_2^2))^2) / (2 * ap^2));

	IGau_Gau_val = M * IGau(bt, ap, Min_Range, Max_Range);

	return IGau_Gau_val;
}

double IGau(double Mu, double StDev, double Min_Range, double Max_Range)
{
	double C, ap, IGau_val;

	C = -1.0 / (2.0 * pow(StDev, 2));

	ap = 1.72254 / pow(StDev, 0.9772);
	IGau_val = sqrt(2 * 3.1415927) * StDev * (1.0 / (1.0 + exp(-1.0 * ap * (Max_Range - Mu))) - 1.0 / (1.0 + exp(-1.0 * ap * (Min_Range - Mu))));
//	IGau_val = sqrt(2 * 3.1415927) * StDev * (1.0 / (1.0 + Exponential(ap * (Max_Range - Mu))) - 1.0 / (1.0 + Exponential(ap * (Mu - Min_Range))));

	return IGau_val;
}

int Matrix_Inversion(double **Jacobian, double **Inv_Jacobian, int order, double *det)
{
	int npivot;		//No. of times rows are interchange
	int pass, row, col, maxrow, i, j, error_flag;
	double temp, pivot, mult;

	//Initialization
	//Store the identity matrix in Inv_Jacobian

	for (i=0; i< order; ++i)
	{
		for(j=0; j< order; ++j)
		{
			if (i == j)
				Inv_Jacobian[i][j] = 1.0;
			else
				Inv_Jacobian[i][j] = 0.0;
		}
	}

	*det = 1.0;
	npivot = 0;

	for (pass = 0; pass < order; ++pass)
	{
		maxrow = pass;
		for (row = pass; row < order; ++row)
			if (fabs(Jacobian[row][pass]) > fabs(Jacobian[maxrow][pass]))
				maxrow = row;

		if (maxrow != pass)
			++npivot;

		for (col = 0; col < order; ++col)
		{
			temp = Jacobian[pass][col];
			Jacobian[pass][col] = Jacobian[maxrow][col];
			Jacobian[maxrow][col] = temp;

			temp = Inv_Jacobian[pass][col];
			Inv_Jacobian[pass][col] = Inv_Jacobian[maxrow][col];
			Inv_Jacobian[maxrow][col] = temp;
		}

		pivot = Jacobian[pass][pass];
		*det *= pivot;

		if (fabs(*det) < 1.0e-40)
		{
			error_flag = 0;
			return(error_flag);
		}

		for(col = 0; col < order; ++col)
		{
			Jacobian[pass][col] = Jacobian[pass][col] / pivot;
			Inv_Jacobian[pass][col] = Inv_Jacobian[pass][col] / pivot;
		}

		for(row = 0; row < order; ++row)
		{
			if (row != pass)
			{
				mult = Jacobian[row][pass];
				for (col = 0; col < order; ++col)
				{
					Jacobian[row][col] = Jacobian[row][col] - Jacobian[pass][col] * mult;
					Inv_Jacobian[row][col] = Inv_Jacobian[row][col] - Inv_Jacobian[pass][col] * mult;
				}
			}
		}
	}

	if (npivot % 2 != 0)
		*det = *det * -1.0;

	return 1;
}

void OrthogonalMatrix_Generator(int no_dimension, double ***O_Matrix)
{
	int i, j, element_index;
	double *a, *x, x_mag, prd_cos;

	a = (double *)malloc(sizeof(double) * no_dimension);
	x = (double *)malloc(sizeof(double) * no_dimension);

	for(i=0; i< no_dimension-1; ++i)
		a[i] = (double) rand() / RAND_MAX * 2.0 * PI;
	a[no_dimension-1] = 0.0;

	x[0] = sin(a[0]);
	x_mag = pow(x[0],2);
	prd_cos = cos(a[0]);
	for(i=1; i< no_dimension; ++i)
	{
		x[i] = prd_cos * sin(a[i]);
		x_mag += pow(x[i], 2);
		prd_cos *= cos(a[i]);
	}

	for(i=0; i< no_dimension; ++i)
		for(j=0; j< no_dimension; ++j)
			(*O_Matrix)[i][j] = -2.0 * x[i] * x[j] / x_mag;

	for(i=0; i< no_dimension; ++i)
		(*O_Matrix)[i][i] += 1.0;

	free(a);
	free(x);

	for(i=0; i< no_dimension; ++i)
	{
		element_index = rand() % no_dimension;

		for(j=0; j< no_dimension; ++j)
			(*O_Matrix)[element_index][j] *= -1.0;

		for(j=0; j< no_dimension; ++j)
			(*O_Matrix)[j][element_index] *= -1.0;
	}

	return;
}

void Orthonormal_Basis(double *Base, int no_Dimension, double ***OrthonormalBasis)
{
	int i, j, k;
	double *cur_Base, magnitude, weight;

	cur_Base = (double *)malloc(sizeof(double) * no_Dimension);

	magnitude = 0.0;
	for(k=0; k< no_Dimension; ++k)
	{
		cur_Base[k] = Base[k];
		magnitude += pow(cur_Base[k], 2);
	}
	magnitude = sqrt(magnitude);

	for(k=0; k< no_Dimension; ++k)
		cur_Base[k] /= magnitude;

	for(i =0; i < no_Dimension - 1; ++i)
	{
		for(k=0; k< no_Dimension; ++k)
			(*OrthonormalBasis)[i][k] = (double) rand() / RAND_MAX;

		weight = 0.0;
		for(k=0; k< no_Dimension; ++k)
			weight += cur_Base[k] * (*OrthonormalBasis)[i][k];
		for(k=0; k< no_Dimension; ++k)
			(*OrthonormalBasis)[i][k] -= weight * cur_Base[k];

		for(j =0; j < i; ++j)
		{
			weight = 0.0;
			for(k=0; k< no_Dimension; ++k)
				weight += (*OrthonormalBasis)[j][k] * (*OrthonormalBasis)[i][k];
			for(k=0; k< no_Dimension; ++k)
				(*OrthonormalBasis)[i][k] -= weight * (*OrthonormalBasis)[j][k];
		}

		magnitude = 0.0;
		for(k=0; k< no_Dimension; ++k)
			magnitude += pow((*OrthonormalBasis)[i][k], 2);
		magnitude = sqrt(magnitude);

		for(k=0; k< no_Dimension; ++k)
			(*OrthonormalBasis)[i][k] /= magnitude;
	}

	free(cur_Base);

	return;
}

double Min(double *data_list, int no_data)
{
	int i;
	double min_data;

	min_data = data_list[0];
	for(i=1; i< no_data; ++i)
		if(min_data > data_list[i])
			min_data = data_list[i];

	return min_data;
}

double Max(double *data_list, int no_data)
{
	int i;
	double max_data;

	max_data = data_list[0];
	for(i=1; i< no_data; ++i)
		if(max_data < data_list[i])
			max_data = data_list[i];

	return max_data;
}

double Sigma(double x, double ap)
{
	double sigma;

	sigma = 2.0 / (1.0 + exp(-1.0 * ap * x)) - 1.0;

	return sigma;
}

double dSigma(double x, double ap)
{
	double exp_x, dsigma;

	exp_x = exp(-1.0 * ap * x);

//	dsigma = 2.0 * ap * exp_x / pow(1.0 + exp_x, 2.0);
	dsigma = 2.0 * ap / (exp_x * pow(1.0 / exp_x + 1.0, 2.0));

	return dsigma;
}


double nChoosek(int n, int k)
{
	int i;
	double value1, value2;

	value1 = value2 = 1;

	for(i=k+1; i<= n; ++i)
		value1 *= (double) i;

	for(i=1; i<= n-k; ++i)
		value2 *= (double) i;

	return value1 / value2;
}

void Dec2Bin(int dec, double **bin, int no_digit)
{
	int i;

	for(i=0; i< no_digit; ++i)
	{
		(*bin)[i] = dec & 1;
		dec >>= 1;
	}

	return;
}

void Pause()
{
	int i;

	printf("Pause ......");
	scanf("%d", &i);

	return;
}

int Sign(double val)
{
	if(val >= 0.0)
		return 1;

	return -1;
}

