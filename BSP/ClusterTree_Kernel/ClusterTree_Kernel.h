#ifndef CLUSTERTREE_KERNEL

#define CLUSTERTREE_KERNEL
#include <stdio.h>
// * = 1-D array without fixed dimension
// ** = 2-D array without fixed dimension

struct CLUSTER_NODE {
	int ID_Self;			//The ID of the Current Node
	int ID_Parent;			//The ID of the Parent Node

	double Optimal_Fitness;		 //The optimal fitness in history
	double *Optimal_Data;		 //The optimal data in history (the current node information)

	double no_Visit;		//No. of Visit
	int Dimension;			//Dimension Size
	// v2
	int pseudo_Dimension;			//Dimension Size
	double Threshold;		// <= ?

	double Survival_Rate;

	struct CLUSTER_NODE *Parent;
	struct CLUSTER_NODE *Child_Left;
	struct CLUSTER_NODE *Child_Right;

	//------------ For 'Blockader'
	int Blockader;

	//------------ For 'Local Basin'
	struct CLUSTER_NODE *Basin_Node;
};

typedef struct CLUSTER_NODE Cluster_Node;

typedef struct CLUSTER_TREE {
	int Current_ID;

	int no_Dimension;
	int no_Leaf;			//number of node inside the tree
	double **Interval; 		//array storing the whole search space

	int max_Archive_Size;
	int no_AxisBoundary;
	int *AxisBoundary_Set;

	Cluster_Node *Root;
	double **cur_Interval;	//array storing the interval of current focusing node
} Cluster_Tree;


void Cluster_Tree_Insertion(Cluster_Tree*, double*, double, Cluster_Node**);
void Cluster_Tree_SearchInterval_Set(Cluster_Tree*);
void Cluster_Tree_Write(Cluster_Tree*, char*);
void Cluster_Tree_Construction(Cluster_Tree*, int, int, int);
void Cluster_Tree_Destruction(Cluster_Tree);

void Cluster_Node_Path(Cluster_Node*, int*, int);
void Cluster_Node_Insertion(Cluster_Node*, double**, double*, double, int*, int, int*, int**, int*);
void Cluster_Node_Interval(Cluster_Node*, Cluster_Tree*);
void Cluster_Node_Write(Cluster_Node*, double***, int, FILE*);
void Cluster_Node_Contruction(Cluster_Node*, Cluster_Node*, int);
void Cluster_Node_Destruction(Cluster_Node);

double Cluster_Node_Difference(Cluster_Node*, double*, int);

Cluster_Node *Cluster_Node_Search(Cluster_Node*, double***, double*);

// v2
Cluster_Node *Cluster_Node_Search_for_Insert(Cluster_Node *root_Node, double ***cur_Interval, double *cur_Data);




#endif
