#pragma once
#include <vector>
#include "../ClusterTree_Kernel/ClusterTree_Kernel.h"
// Roulette
struct Selection_Strategy
{
	Cluster_Node *cur_Node;
	double probability;
	double accum_probability;
};
void Archive_Leaf(Cluster_Node *cur_Node, std::vector<Selection_Strategy> &Select_Vector);
void Cal_Archive_Leaf_Probility(std::vector<Selection_Strategy> &Leaf_Vector);
Cluster_Node* Selection_Roulette(const double &rand, std::vector<Selection_Strategy> &Leaf_Vector);
