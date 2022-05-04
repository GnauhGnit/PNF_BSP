/*
* Archive.h
*
*  Created on: Jul 13, 2017
*      Author: hting
*/

#ifndef ARCHIVE_H_
#define ARCHIVE_H_
#include <vector>
#include "../header.h"
void Archive(double *individual, double fitness, std::vector<std::vector<double> > &Archive,
	std::vector<double> &fitnessArchive, double &maxFitness, int dim, double radius);
void Final_Archive_4_Evaluation_V3(std::vector<std::vector<double> > &Archive,
	std::vector<double> &fitnessArchive, std::vector<std::vector<double> > &Final_Archive, std::vector<double> &Final_fitnessArchive,
	struct NewType *fitness_index,
	const double &maxFitness, const double &epsilon, const double &radius, const int &dim);
#endif /* ARCHIVE_H_ */
