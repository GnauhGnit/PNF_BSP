/******************************************************************************************/
//	Program Name	:	Non-Revisiting Genetic Algorithm (NrGA)
//	File Name		:	ClusterTree_Kernel.c
//	Auther			:	Dr. Chow Chi Kin, Dr. Yuen Shiu Yin
//	Edit by			:	Leung Shing Wa
//	University		:	City University of Hong Kong
//	Department		:	Electronic Engineering
//	Last Update		:	10 Sep 2009
//	Reference		:	A Genetic Algorithm that Adaptively Mutates and Never Revisits, 
//						IEEE Transactions on Evolutionary Computation, 
//						Vol 13(2) (April 2009) 454-472.
//
//	Discription		:	Basic functions for BSP Tree are defined in this file.
#include "../../BSP/ClusterTree_Kernel/ClusterTree_Kernel.h"

#include<stdio.h>
#include<stdlib.h>
#include<malloc.h>
#include<math.h>

//Function for inserting data (chromosomes) into the Tree
void Cluster_Tree_Insertion(Cluster_Tree *cur_Tree, double *cur_Data, double cur_Fitness, Cluster_Node **target_Node)
{
	double node_difference;
	Cluster_Node *cur_Node;

	// Classifer the current data belongs which cluster, 2 things will be returned
	// cur_Node = the cluster containing current data
	// cur_Interval = the interval of the cluster
	// 如果只有根节点，直接返回。如果根节点还有孩子节点，则继续搜索
//	cur_Node = Cluster_Node_Search((*cur_Tree).Root, &(*cur_Tree).cur_Interval, cur_Data);   //Interval = lower bound& upper bound
// v2
		cur_Node = Cluster_Node_Search_for_Insert((*cur_Tree).Root, &(*cur_Tree).cur_Interval, cur_Data);   //Interval = lower bound& upper bound

	// 计算cur_Node与cur_Data之间的差异
	// Calculate the distance between the node and the data
	node_difference = Cluster_Node_Difference(cur_Node, cur_Data, (*cur_Tree).no_Dimension);

	// if Current Tree have not parent (it is root node) or distance > 0 (it is not revisited)
	if((*cur_Node).Parent == NULL || node_difference > 0)	
	{
		// insert the current data to the tree
		Cluster_Node_Insertion(cur_Node, (*cur_Tree).cur_Interval, cur_Data, cur_Fitness, &(*cur_Tree).no_Leaf, (*cur_Tree).no_Dimension, &(*cur_Tree).no_AxisBoundary, &(*cur_Tree).AxisBoundary_Set, &(*cur_Tree).Current_ID);

		if((*cur_Node).Parent != NULL)
		{
			*target_Node = cur_Node;

			return;
		}
	}

	// for the case revisiting occur
	*target_Node = NULL;

	return;
}

//Function to read the target interval set to the Tree
void Cluster_Tree_SearchInterval_Set(Cluster_Tree *cur_Tree)
{
	int i;

	for(i=0; i< (*cur_Tree).no_Dimension; ++i)
	{
		(*cur_Tree).cur_Interval[i][0] = (*cur_Tree).Interval[i][0];
		(*cur_Tree).cur_Interval[i][1] = (*cur_Tree).Interval[i][1];
	}

	return;
}

//Function to print out the Tree
void Cluster_Tree_Write(Cluster_Tree *cur_Tree, char *tree_filename)
{
	FILE *tree_ptr;

	tree_ptr = fopen(tree_filename, "w");

	fprintf(tree_ptr, "%d\n", (*cur_Tree).no_Dimension);

	Cluster_Tree_SearchInterval_Set(cur_Tree);

	Cluster_Node_Write((*cur_Tree).Root, &(*cur_Tree).cur_Interval, (*cur_Tree).no_Dimension, tree_ptr);

	fclose(tree_ptr);

	return;
}

// Function to initialize the Tree
void Cluster_Tree_Construction(Cluster_Tree *cur_Tree, int cur_Dimension, int interval_Size, int max_archive_size)
{
	int i;

	(*cur_Tree).no_Dimension = cur_Dimension;
	(*cur_Tree).no_Leaf = 0;
	(*cur_Tree).max_Archive_Size = max_archive_size;
	(*cur_Tree).Root = (Cluster_Node *)malloc(sizeof(Cluster_Node));
	(*(*cur_Tree).Root).Dimension = -1;
	(*(*cur_Tree).Root).pseudo_Dimension = -1;
	(*(*cur_Tree).Root).Parent = NULL;
	(*(*cur_Tree).Root).Optimal_Data = (double *)malloc(sizeof(double) * cur_Dimension);
	(*(*cur_Tree).Root).Child_Left = NULL;
	(*(*cur_Tree).Root).Child_Right = NULL;
	(*(*cur_Tree).Root).no_Visit = 0;
	(*cur_Tree).Current_ID = 0;

	(*cur_Tree).AxisBoundary_Set = (int *)malloc(sizeof(int) * cur_Dimension);
	(*cur_Tree).Interval = (double **)malloc(sizeof(double*) * cur_Dimension);
	(*cur_Tree).cur_Interval = (double **)malloc(sizeof(double*) * cur_Dimension);
	for(i=0; i< cur_Dimension; ++i)
	{
		(*cur_Tree).Interval[i] = (double *)malloc(sizeof(double) * 2);
		(*cur_Tree).cur_Interval[i] = (double *)malloc(sizeof(double) * 2);
	}

	return;
}
void Cluster_Node_Destruction(Cluster_Node *cur_Node)
{
	if(cur_Node == NULL)
	{
		return;
	}
	Cluster_Node_Destruction((*cur_Node).Child_Left);
	Cluster_Node_Destruction((*cur_Node).Child_Right);
	free(cur_Node);
}
// Function to delete the tree
void Cluster_Tree_Destruction(Cluster_Tree cur_Tree)
{
	int i;

	Cluster_Node_Destruction(*cur_Tree.Root);

	free(cur_Tree.Root);


	for(i=0; i< cur_Tree.no_Dimension; ++i)
		free(cur_Tree.cur_Interval[i]);
	free(cur_Tree.cur_Interval);

	return;

}

// Function to search the cluster which containing the current data
Cluster_Node *Cluster_Node_Search(Cluster_Node *root_Node, double ***cur_Interval, double *cur_Data)
{
	Cluster_Node *cur_Node;

	cur_Node = root_Node;
	while ((*cur_Node).Dimension > -1)
	{
		if(cur_Data[(*cur_Node).Dimension] < (*cur_Node).Threshold)
		{
			(*cur_Interval)[(*cur_Node).Dimension][1] = (*cur_Node).Threshold;
			cur_Node = (*cur_Node).Child_Left;
		}
		else
		{
			(*cur_Interval)[(*cur_Node).Dimension][0] = (*cur_Node).Threshold;
			cur_Node = (*cur_Node).Child_Right;
		}
	}

	return cur_Node;
}
// v2
Cluster_Node *Cluster_Node_Search_for_Insert(Cluster_Node *root_Node, double ***cur_Interval, double *cur_Data)
{
	Cluster_Node *cur_Node;

	cur_Node = root_Node;
	while ((*cur_Node).Dimension > -1)
	{
		if(cur_Data[(*cur_Node).Dimension] < (*cur_Node).Threshold)
		{
			(*cur_Interval)[(*cur_Node).Dimension][1] = (*cur_Node).Threshold;
			cur_Node = (*cur_Node).Child_Left;
		}
		else
		{
			(*cur_Interval)[(*cur_Node).Dimension][0] = (*cur_Node).Threshold;
			cur_Node = (*cur_Node).Child_Right;
		}
	}

	return cur_Node;
}
//function to insert a new node into the tree
void Cluster_Node_Insertion(Cluster_Node *leaf_Node, double **cur_Interval, double *cur_Data, double cur_Fitness, int *no_Leaf, int cur_Dimension, int *no_AxisBoundary, int **AxisBoundary_Set, int *Current_ID)
{
	int i;
	double interval_diff;

	++(*no_Leaf);	//increase the number of leaf node

	if(leaf_Node == NULL)
		printf("NULL\n");

	//当前节点为根节点，直接插入替换
	if((*leaf_Node).no_Visit == 0 && (*leaf_Node).Parent == NULL)
	{
		//============ Root Node
		++(*leaf_Node).no_Visit;
		for(i=0; i< cur_Dimension; ++i)
			(*leaf_Node).Optimal_Data[i] = cur_Data[i];
		(*leaf_Node).Optimal_Fitness = cur_Fitness;

		(*leaf_Node).ID_Parent = -1;
		(*leaf_Node).ID_Self = (*Current_ID);
		++(*Current_ID);
	}
	//非根节点
	else
	{	//============ Not a Root Node
		for(i=0; i< cur_Dimension; ++i)
			if((*leaf_Node).Optimal_Data[i] != cur_Data[i] &&
			   ((*leaf_Node).Dimension == -1 || interval_diff < fabs((*leaf_Node).Optimal_Data[i] - cur_Data[i])))
			{
				(*leaf_Node).Dimension = i;
				(*leaf_Node).pseudo_Dimension = i;
				//Assume it is used in integer, Threshold = the critical line
				(*leaf_Node).Threshold = ceil(0.5 * ((*leaf_Node).Optimal_Data[i] + cur_Data[i])) - 0.5;
				interval_diff = fabs((*leaf_Node).Optimal_Data[i] - cur_Data[i]);
			}

		(*leaf_Node).Child_Left = (Cluster_Node *)malloc(sizeof(Cluster_Node));
		(*leaf_Node).Child_Right = (Cluster_Node *)malloc(sizeof(Cluster_Node));

		Cluster_Node_Contruction(leaf_Node, (*leaf_Node).Child_Right, cur_Dimension);
		Cluster_Node_Contruction(leaf_Node, (*leaf_Node).Child_Left, cur_Dimension);

		if((*leaf_Node).Optimal_Data[(*leaf_Node).Dimension] < cur_Data[(*leaf_Node).Dimension])
		{
			for(i=0; i< cur_Dimension; ++i)
			{
				(*(*leaf_Node).Child_Left).Optimal_Data[i] = (*leaf_Node).Optimal_Data[i];
				(*(*leaf_Node).Child_Right).Optimal_Data[i] = cur_Data[i];
			}
			(*(*leaf_Node).Child_Left).Optimal_Fitness = (*leaf_Node).Optimal_Fitness;
			(*(*leaf_Node).Child_Right).Optimal_Fitness = cur_Fitness;

			(*(*leaf_Node).Child_Left).ID_Parent = (*(*leaf_Node).Child_Right).ID_Parent = (*leaf_Node).ID_Self;

			(*(*leaf_Node).Child_Left).ID_Self = (*Current_ID);
			(*(*leaf_Node).Child_Right).ID_Self = (*Current_ID) + 1;
			(*Current_ID) += 2;
		}
		else
		{
			for(i=0; i< cur_Dimension; ++i)
			{
				(*(*leaf_Node).Child_Left).Optimal_Data[i] = cur_Data[i];
				(*(*leaf_Node).Child_Right).Optimal_Data[i] = (*leaf_Node).Optimal_Data[i];
			}
			(*(*leaf_Node).Child_Left).Optimal_Fitness = cur_Fitness;
			(*(*leaf_Node).Child_Right).Optimal_Fitness = (*leaf_Node).Optimal_Fitness;

			(*(*leaf_Node).Child_Left).ID_Parent = (*(*leaf_Node).Child_Right).ID_Parent = (*leaf_Node).ID_Self;

			(*(*leaf_Node).Child_Left).ID_Self = (*Current_ID);
			(*(*leaf_Node).Child_Right).ID_Self = (*Current_ID) + 1;
			(*Current_ID) += 2;
		}
	}

	return;
}

//Function to copy the interval of the target node to cur_interval
void Cluster_Node_Interval(Cluster_Node *tgr_Node, Cluster_Tree *cur_Tree)
{
	Cluster_Node *cur_Node;

	Cluster_Tree_SearchInterval_Set(cur_Tree);

	cur_Node = (*cur_Tree).Root;
	while (cur_Node != tgr_Node && (*cur_Node).Dimension > -1)
	{
		if((*tgr_Node).Optimal_Data[(*cur_Node).Dimension] < (*cur_Node).Threshold)
		{
			(*cur_Tree).cur_Interval[(*cur_Node).Dimension][1] = (*cur_Node).Threshold;
			cur_Node = (*cur_Node).Child_Left;
		}
		else
		{
			(*cur_Tree).cur_Interval[(*cur_Node).Dimension][0] = (*cur_Node).Threshold;
			cur_Node = (*cur_Node).Child_Right;
		}
	}

	return;
}

//Function for update the tree.
void Cluster_Node_Path(Cluster_Node *leaf_Node, int *no_Leaf, int cur_Dimension)
{
	int i;
	Cluster_Node *parent_Node, *child_Node;

	child_Node = leaf_Node;
	parent_Node = (*leaf_Node).Parent;
	while (parent_Node != NULL)	// it will stop when the pointer is pointing to the root node
	{
		++(*parent_Node).no_Visit;

		//update the optimal data find in the area under this node
		if((*parent_Node).Optimal_Fitness < (*child_Node).Optimal_Fitness)
		{
			(*parent_Node).Optimal_Fitness = (*child_Node).Optimal_Fitness;
			for(i=0; i< cur_Dimension; ++i)
				(*parent_Node).Optimal_Data[i] = (*child_Node).Optimal_Data[i];
		}

		//prune the tree if all the area under this node is visited 
		if((*(*parent_Node).Child_Left).Dimension == -2 && (*(*parent_Node).Child_Right).Dimension == -2)
		{
			--(*no_Leaf);

			free((*(*parent_Node).Child_Left).Optimal_Data);
			free((*parent_Node).Child_Left);
			free((*(*parent_Node).Child_Right).Optimal_Data);
			free((*parent_Node).Child_Right);

			(*parent_Node).Dimension = -2;
			(*parent_Node).Child_Left = NULL;
			(*parent_Node).Child_Right = NULL;
		}

		child_Node = parent_Node;
		parent_Node = (*parent_Node).Parent;
	}

	return;
}

//function to calculate the distance between the data and the node (to check revist)
double Cluster_Node_Difference(Cluster_Node *cur_Node, double *cur_Data, int no_dimension)
{
	int i;
	double difference;

	difference = 0.0;
	for(i=0; i< no_dimension; ++i)
		difference += pow((*cur_Node).Optimal_Data[i] - cur_Data[i], 2);

	return sqrt(difference);
}

void Cluster_Node_Write(Cluster_Node *cur_Node, double ***cur_Interval, int no_dimension, FILE *tree_ptr)
{
	int i;
	double min_bound, max_bound;

	if((*cur_Node).Dimension > -1)
	{
		min_bound = (*cur_Interval)[(*cur_Node).Dimension][0];
		max_bound = (*cur_Interval)[(*cur_Node).Dimension][1];

		(*cur_Interval)[(*cur_Node).Dimension][0] = min_bound;
		(*cur_Interval)[(*cur_Node).Dimension][1] = (*cur_Node).Threshold;
		Cluster_Node_Write((*cur_Node).Child_Left, cur_Interval, no_dimension, tree_ptr);

		(*cur_Interval)[(*cur_Node).Dimension][0] = (*cur_Node).Threshold;
		(*cur_Interval)[(*cur_Node).Dimension][1] = max_bound;
		Cluster_Node_Write((*cur_Node).Child_Right, cur_Interval, no_dimension, tree_ptr);
		
		(*cur_Interval)[(*cur_Node).Dimension][0] = min_bound;
		(*cur_Interval)[(*cur_Node).Dimension][1] = max_bound;

		for(i=0; i< no_dimension; ++i)
			fprintf(tree_ptr, "%lf %lf ", (*cur_Interval)[i][0], (*cur_Interval)[i][1]);
//		fprintf(tree_ptr, "%1.0lf %d %d\n", (*cur_Node).no_Visit, (*cur_Node).ID_Self, (*cur_Node).ID_Parent);
		fprintf(tree_ptr, "%1.0lf %d %d (%d, %d)\n", (*cur_Node).no_Visit, (*cur_Node).ID_Self, (*cur_Node).ID_Parent, (*cur_Node).Dimension, (*cur_Node).pseudo_Dimension);
	}
	else
	{
		for(i=0; i< no_dimension; ++i)
			fprintf(tree_ptr, "%lf %lf ", (*cur_Interval)[i][0], (*cur_Interval)[i][1]);
//		fprintf(tree_ptr, "%1.0lf %d %d\n", (*cur_Node).no_Visit, (*cur_Node).ID_Self, (*cur_Node).ID_Parent);
		fprintf(tree_ptr, "%1.0lf %d %d (%d, %d)\n", (*cur_Node).no_Visit, (*cur_Node).ID_Self, (*cur_Node).ID_Parent, (*cur_Node).Dimension, (*cur_Node).pseudo_Dimension);
	}

	return;
}

//function to build and initializate the cluster node
void Cluster_Node_Contruction(Cluster_Node *parent_Node, Cluster_Node *child_Node, int cur_Dimension)
{
	(*child_Node).ID_Self = -1;
	(*child_Node).ID_Parent = -1;
	(*child_Node).Dimension = -1;			//Dimension  = 1 means it is a void node
	// v2
	(*child_Node).pseudo_Dimension = -1;			//Dimension  = 1 means it is a void nodes
	(*child_Node).Parent = parent_Node;
	(*child_Node).Child_Left = NULL;
	(*child_Node).Child_Right = NULL;
	(*child_Node).no_Visit = 1;
	(*child_Node).Optimal_Data = (double *)malloc(sizeof(double) * cur_Dimension);
	(*child_Node).Survival_Rate = 0;
	(*child_Node).Blockader = -1;

	return;
}

//function to destruction the cluster node
void Cluster_Node_Destruction(Cluster_Node cur_Node)
{

	if(cur_Node.Child_Left != NULL)
	{
		Cluster_Node_Destruction(*cur_Node.Child_Left);
		free(cur_Node.Child_Left);
	}

	if(cur_Node.Child_Right != NULL)
	{
		Cluster_Node_Destruction(*cur_Node.Child_Right);
		free(cur_Node.Child_Right);
	}

	free(cur_Node.Optimal_Data);

	return;
}
