#include "./Selection_Strategy.h"
#include "../../header.h"
using namespace std;


void Archive_Leaf(Cluster_Node *cur_Node, vector<Selection_Strategy> &Select_Vector)
{
	Selection_Strategy node;
	node.cur_Node = cur_Node;
	node.probability = 0.0;
	node.accum_probability = 0.0;
	Select_Vector.push_back(node);
}

void Cal_Archive_Leaf_Probility(vector<Selection_Strategy> &Leaf_Vector)
{
	double fitness = 0.0;
	int no_visit = 0;
	double total_prob = 0.0;
	double acum_prob = 0.0;
	double *Fitness_Leaf = new double[Leaf_Vector.size()];
	double *Visit_Leaf = new double[Leaf_Vector.size()];
	double min_Fitness = 0.0;
	int min_visit = 0;
	for (size_t i = 0; i < Leaf_Vector.size(); i++)
	{
		Fitness_Leaf[i] = (*Leaf_Vector[i].cur_Node).Optimal_Fitness;
		Visit_Leaf[i] = (*Leaf_Vector[i].cur_Node).no_Visit;
		min_Fitness = ((i == 0) || Fitness_Leaf[i] < min_Fitness) ? Fitness_Leaf[i] : min_Fitness;
		min_visit = ((i == 0) || Visit_Leaf[i] < min_visit) ? Visit_Leaf[i] : min_visit;
	}
	for (size_t i = 0; i < Leaf_Vector.size(); i++)
	{
		Leaf_Vector[i].probability = (Fitness_Leaf[i] - min_Fitness) / (double)(Visit_Leaf[i] - min_visit + 1);
		total_prob += Leaf_Vector[i].probability;
	}
	for (size_t i = 0; i < Leaf_Vector.size(); i++)
	{
		Leaf_Vector[i].probability /= total_prob;
		acum_prob += Leaf_Vector[i].probability;
		Leaf_Vector[i].accum_probability = acum_prob;
	}
	//for (size_t i = 0; i < Leaf_Vector.size(); i++)
	//{
	//	fitness = (*Leaf_Vector[i].cur_Node).Optimal_Fitness;
	//	no_visit = (*Leaf_Vector[i].cur_Node).no_Visit;
	//	Leaf_Vector[i].probability = fitness / (double)no_visit;
	//	total_prob += Leaf_Vector[i].probability;
	//}
	//for (size_t i = 0; i < Leaf_Vector.size(); i++)
	//{
	//	Leaf_Vector[i].probability /= total_prob;
	//	acum_prob += Leaf_Vector[i].probability;
	//	Leaf_Vector[i].accum_probability = acum_prob;
	//}
	delete[]Fitness_Leaf;
	delete[]Visit_Leaf;
}
Cluster_Node* Selection_Roulette(const double &rand, vector<Selection_Strategy> &Leaf_Vector)
{
	for (size_t i = 0; i < Leaf_Vector.size(); i++)
	{
		if (rand < Leaf_Vector[i].accum_probability)
			return Leaf_Vector[i].cur_Node;
	}
	return NULL;
}
