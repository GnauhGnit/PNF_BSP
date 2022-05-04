#ifndef NONREVISITINGSCHEME_KERNEL

#define NONREVISITINGSCHEME_KERNEL

#ifndef CLUSTERTREE_KERNEL
	#include"../../BSP/ClusterTree_Kernel/ClusterTree_Kernel.h"
#endif


typedef struct NONREVISITINGSCHEME {
	int Activate;

	int Max_no_Leaf;

	Cluster_Tree Solution_Tree;

//	AutoPrunning Prunning_Info;		
} NonRevisitingScheme;


//============ UnConstraint Optimization
void NonRevisitingScheme_TreePruning(Cluster_Node *cur_Node, int *no_Leaf, int cur_Dimension);
void NonRevisitingScheme_UnConstraint_Update(NonRevisitingScheme*, double*, double);
int  NonRevisitingScheme_UnConstraint_SaturateVerification(double**, int);
void NonRevisitingScheme_UnConstraint_TreePath(Cluster_Node*, double**, int*, int);
void NonRevisitingScheme_UnConstraint_NodeSearch(double**, int, Cluster_Tree*, Cluster_Node**);

#endif