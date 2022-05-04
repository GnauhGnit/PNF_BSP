/******************************************************************************************/
//	Program Name	:	Non-Revisiting Genetic Algorithm (NrGA)
//	File Name		:	NonRevisitingScheme_Kernel.c
//	Auther			:	Dr. Chow Chi Kin, Dr. Yuen Shiu Yin
//	Edit by			:	Leung Shing Wa
//	University		:	City University of Hong Kong
//	Department		:	Electronic Engineering
//	Last Update		:	10 Sep 2009
//	Reference		:	A Genetic Algorithm that Adaptively Mutates and Never Revisits, 
//						IEEE Transactions on Evolutionary Computation, 
//						Vol 13(2) (April 2009) 454-472.
//
//	Discription		:	The functions to maintain (update / use) the BSP tree can be found in this file.
#include "../../BSP/NonRevisitingScheme_Kernel/NonRevisitingScheme_Kernel.h"

#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#include<malloc.h>
#include<math.h>
#include<time.h>

#include "../../BSP/ACC_Math_Kernel/ACC_Math_Kernel.h"
#include "../ClusterTree_Kernel/ClusterTree_Kernel.h"


/*The TreePruning fucntion, it will prune the space that is fully filled*/
void NonRevisitingScheme_TreePruning(Cluster_Node *cur_Node, int *no_Leaf, int cur_Dimension)
{// 从叶子节点往上找，如果左右节点都为-2，则父节点也设置为-2.即没有剩余搜索空间

	int i;

	while (cur_Node != NULL)
	{
		++(*cur_Node).no_Visit;

//		//用左右节点更新父节点的Optimal_Fitness
//		if((*cur_Node).Optimal_Fitness > (*(*cur_Node).Child_Left).Optimal_Fitness)
//		{
//			(*cur_Node).Optimal_Fitness = (*(*cur_Node).Child_Left).Optimal_Fitness;
//			for(i=0; i< cur_Dimension; ++i)
//				(*cur_Node).Optimal_Data[i] = (*(*cur_Node).Child_Left).Optimal_Data[i];
//		}
//		else
//		{
//			(*cur_Node).Optimal_Fitness = (*(*cur_Node).Child_Right).Optimal_Fitness;
//			for(i=0; i< cur_Dimension; ++i)
//				(*cur_Node).Optimal_Data[i] = (*(*cur_Node).Child_Right).Optimal_Data[i];
//		}

		// v2
		if((*(*cur_Node).Child_Left).Dimension == -2 && (*(*cur_Node).Child_Right).Dimension == -2)
		{
			(*cur_Node).Dimension = -2;			// Dimension  = -2 means it is leaf node
		}

//		if((*(*cur_Node).Child_Left).Dimension == -2 && (*(*cur_Node).Child_Right).Dimension == -2)
//		{
//			--(*no_Leaf);
//
//			free((*(*cur_Node).Child_Left).Optimal_Data);
//			free((*cur_Node).Child_Left);
//			free((*(*cur_Node).Child_Right).Optimal_Data);
//			free((*cur_Node).Child_Right);
//
//			(*cur_Node).Dimension = -2;			// Dimension  = -2 means it is leaf node
//			(*cur_Node).Child_Left = NULL;
//			(*cur_Node).Child_Right = NULL;
//		}

		cur_Node = (*cur_Node).Parent;
	}

	return;
}

/****************************
 *
 *	Function: Update the cluster tree. It will insert the data into the cluster tree by calling cluster_Tree_Insertion
 *  		and maintain the tree path by calling NonRevisitingScheme Unconstraint Treepath(The treepath function will prune the tree).
 *	Input:  cur_NRS 			(The whole non-revisiting scheme)
 *			cur_Solution      	(The data)
 *			cur_Goodness  		(Fitness of the data)
 *
 */
void NonRevisitingScheme_UnConstraint_Update(NonRevisitingScheme *cur_NRS, double *cur_Solution, double cur_Goodness)
{
	Cluster_Node *target_Node;

	if((*cur_NRS).Activate == 1 && (*(*cur_NRS).Solution_Tree.Root).Dimension != -2)
	{

		// 如果插入的节点为根节点或者是已经访问过的节点，则target为NULL；否则，target为插入的节点
		Cluster_Tree_Insertion(&(*cur_NRS).Solution_Tree, cur_Solution, cur_Goodness, &target_Node);

		// 剪枝检查
		// 一但插入新节点时，新节点如果为NULL,则不进行任何操作直接返回；否则，则检查是否剪枝
		NonRevisitingScheme_UnConstraint_TreePath(target_Node, (*cur_NRS).Solution_Tree.cur_Interval, &(*cur_NRS).Solution_Tree.no_Leaf, (*cur_NRS).Solution_Tree.no_Dimension);

		if((*cur_NRS).Max_no_Leaf < (*cur_NRS).Solution_Tree.no_Leaf)
			(*cur_NRS).Max_no_Leaf = (*cur_NRS).Solution_Tree.no_Leaf;
	}

	return;
}
/************************************
 * 	Function: Search which node in the cluster tree that the input chromosome belongs to.
 * 			  The result will be saved into cur_Node.
 *
 * 	Input: 	  cur_Chromosome           (The input chromosome)
 *			  no_gene                  (the number of genes in each chromosome)
 *			  population_Tree          (The cluster tree)
 *			  cur_Node                 (The result node)
 *
 */
void NonRevisitingScheme_UnConstraint_NodeSearch(double **cur_Chromosome, int no_gene, Cluster_Tree *population_Tree, Cluster_Node **cur_Node)
{
	double min_bound, max_bound;

	//============ Initialization
	// intialize current interval
	Cluster_Tree_SearchInterval_Set(population_Tree);
	
	// 从根节点开始查找
	(*cur_Node) = (*population_Tree).Root;
	do {
		//============ Search the leaf
		while ((*(*cur_Node)).Dimension > -1)
		{// dimension > -1 means that there is leaf node.
			min_bound = (*population_Tree).cur_Interval[(*(*cur_Node)).Dimension][0];
			max_bound = (*population_Tree).cur_Interval[(*(*cur_Node)).Dimension][1];

			if((*cur_Chromosome)[(*(*cur_Node)).Dimension] < (*(*cur_Node)).Threshold)
			{//visit the left node
				(*population_Tree).cur_Interval[(*(*cur_Node)).Dimension][1] = (*(*cur_Node)).Threshold;
				(*cur_Node) = (*(*cur_Node)).Child_Left;
			}
			else
			{// visit the right node
				(*population_Tree).cur_Interval[(*(*cur_Node)).Dimension][0] = (*(*cur_Node)).Threshold;
				(*cur_Node) = (*(*cur_Node)).Child_Right;
			}
		}

		if((*(*cur_Node)).Dimension == -2)
		{// it has been visited.
			//============ Mutate Gene
			if((*cur_Node) == (*(*(*cur_Node)).Parent).Child_Left)
			{// if the node resides in left node, then visits the right node
				(*cur_Chromosome)[(*(*(*cur_Node)).Parent).Dimension] = (*(*(*cur_Node)).Parent).Threshold + 0.5;

				(*population_Tree).cur_Interval[(*(*(*cur_Node)).Parent).Dimension][0] = (*(*(*cur_Node)).Parent).Threshold;
				(*population_Tree).cur_Interval[(*(*(*cur_Node)).Parent).Dimension][1] = max_bound;
				(*cur_Node) = (*(*(*cur_Node)).Parent).Child_Right;
			}
			else
			{// otherwise
				(*cur_Chromosome)[(*(*(*cur_Node)).Parent).Dimension] = (*(*(*cur_Node)).Parent).Threshold - 0.5;

				(*population_Tree).cur_Interval[(*(*(*cur_Node)).Parent).Dimension][0] = min_bound;
				(*population_Tree).cur_Interval[(*(*(*cur_Node)).Parent).Dimension][1] = (*(*(*cur_Node)).Parent).Threshold;
				(*cur_Node) = (*(*(*cur_Node)).Parent).Child_Left;
			}
		}
	} while ((*(*cur_Node)).Dimension > -1);

	return;
}
/******************************
 *
 * 	Function: Make the path to find the best fitness node.
 * 			  It will be called once new node is added.
 * 	Input: 	  leaf_Node          (The input node)
 * 			  cur_Interval       (The input interval)
 * 			  no_Leaf            (The number of leaf nodes in the cluster tree, only for record)
 *			  cur_Dimension 	 (The number of dimension of the GA)
 *
 */
void NonRevisitingScheme_UnConstraint_TreePath(Cluster_Node *leaf_Node, double **cur_Interval, int *no_Leaf, int cur_Dimension)
{
	double min_range, max_range;
	Cluster_Node *parent_Node;

	if(leaf_Node == NULL)
		return;

	parent_Node = leaf_Node;

	min_range = cur_Interval[(*parent_Node).Dimension][0];
	max_range = cur_Interval[(*parent_Node).Dimension][1];

	//============ Verify if the left child is saturated
	cur_Interval[(*parent_Node).Dimension][0] = min_range;
	cur_Interval[(*parent_Node).Dimension][1] = (*parent_Node).Threshold;
	// 如果还有剩余空间，则为叶子节点-1；否则，为闭合节点-2
	(*(*parent_Node).Child_Left).Dimension = NonRevisitingScheme_UnConstraint_SaturateVerification(cur_Interval, cur_Dimension);
//	printf("left: [%f][%f] %d, %d\t", cur_Interval[(*parent_Node).Dimension][0], cur_Interval[(*parent_Node).Dimension][1], (*(*parent_Node).Child_Left).pseudo_Dimension, (*(*parent_Node).Child_Left).Dimension);



	//============ Verify if the right child is saturated
	cur_Interval[(*parent_Node).Dimension][0] = (*parent_Node).Threshold;
	cur_Interval[(*parent_Node).Dimension][1] = max_range;
	// 如果还有剩余空间，则为叶子节点-1；否则，为闭合节点-2
	(*(*parent_Node).Child_Right).Dimension = NonRevisitingScheme_UnConstraint_SaturateVerification(cur_Interval, cur_Dimension);
//	printf("right: [%f][%f] %d, %d\n", cur_Interval[(*parent_Node).Dimension][0], cur_Interval[(*parent_Node).Dimension][1], (*(*parent_Node).Child_Left).pseudo_Dimension, (*(*parent_Node).Child_Left).Dimension);

	//cur_Interval上边的只是暂时用一下。。。。
	cur_Interval[(*parent_Node).Dimension][0] = min_range;
	cur_Interval[(*parent_Node).Dimension][1] = max_range;

		//============ Tree Pruning from 'parent_node'
	// v2
	NonRevisitingScheme_TreePruning(parent_Node, no_Leaf, cur_Dimension);
//	printf("\tparent:  %d, %d\n", (*parent_Node).pseudo_Dimension,  (*parent_Node).Dimension);
	return;
}

//Function to check whether the cur_Interval is full (saturated) or not
int NonRevisitingScheme_UnConstraint_SaturateVerification(double **cur_Interval, int cur_Dimension)
{
	int i;

	for(i=0; i< cur_Dimension; ++i)
		if(cur_Interval[i][1] - cur_Interval[i][0] > 1)
			return -1;

	return -2;
}


