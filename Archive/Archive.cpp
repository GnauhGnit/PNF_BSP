/*
* Archive.cpp
*
*  Created on: Jul 13, 2017
*      Author: hting
*/
#include "./Archive.h"
#include "../header.h"
#include <numeric>
void Archive(double *individual, double fitness, std::vector<std::vector<double> > &Archive,
	std::vector<double> &fitnessArchive, double &minFitness, int dim, double radius)
{
	vector<double> indiv(dim);
	for (int j = 0; j < dim; j++)
		indiv[j] = individual[j];

	Archive.push_back(indiv);
	fitnessArchive.push_back(fitness);
}
// Find the max fitness values for eliminating redundant individuals
void Final_Archive_4_Evaluation_V3(std::vector<std::vector<double> > &Archive,
	std::vector<double> &fitnessArchive, std::vector<std::vector<double> > &Final_Archive, std::vector<double> &Final_fitnessArchive,
	NewType *fitness_index,
	const double &maxFitness, const double &epsilon, const double &radius, const int &dim)
{
	int archive_size = fitnessArchive.size();

	for (int i = 0; i < Final_fitnessArchive.size(); i++)
	{
		Final_Archive[i].clear();
	}
	Final_Archive.clear();
	Final_fitnessArchive.clear();



	bool found;
	double dist;
	int pos;
	for (int i = 0; i < archive_size; i++)
	{
		pos = fitness_index[i].id;
		if (maxFitness - fitnessArchive[pos] > epsilon)
			break;
		Final_Archive.push_back(Archive[pos]);
		Final_fitnessArchive.push_back(fitnessArchive[pos]);

	}

}