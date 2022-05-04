/*
 * RecordRevisitingScheme_Kernel.cpp
 *
 *  Created on: Jul 29, 2017
 *      Author: hting
 */
#include "../../BSP/RecordRevisitingScheme_Kernel/RecordRevisitingScheme_Kernel.h"

#include<stdio.h>
#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#include<malloc.h>
#include<math.h>
#include<time.h>
#include <iostream>
#include <iomanip>
#include "../../header.h"
#include "../SelectionStrategey_Kernel/Selection_Strategy.h"
using namespace std;

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
// v2: 一直找到叶节点
void RecordRevisiting_NodeSearch(double *cur_Chromosome, int no_gene, Cluster_Tree *population_Tree, Cluster_Node **cur_Node)
{
	double min_bound, max_bound;
	//============ Initialization
	Cluster_Tree_SearchInterval_Set(population_Tree);

	// 从根节点开始查找
	(*cur_Node) = NULL;
	(*cur_Node) = (*population_Tree).Root;
	while ((*(*cur_Node)).pseudo_Dimension != -1)
	{
		min_bound = (*population_Tree).cur_Interval[(*(*cur_Node)).pseudo_Dimension][0];
		max_bound = (*population_Tree).cur_Interval[(*(*cur_Node)).pseudo_Dimension][1];

//		if((*cur_Chromosome)[(*(*cur_Node)).Dimension] < (*(*cur_Node)).Threshold)
		if(cur_Chromosome[(*(*cur_Node)).pseudo_Dimension] < (*(*cur_Node)).Threshold)
		{//visit the left node
			(*population_Tree).cur_Interval[(*(*cur_Node)).pseudo_Dimension][1] = (*(*cur_Node)).Threshold;
			(*cur_Node) = (*(*cur_Node)).Child_Left;
		}
		else
		{// visit the right node
			(*population_Tree).cur_Interval[(*(*cur_Node)).pseudo_Dimension][0] = (*(*cur_Node)).Threshold;
			(*cur_Node) = (*(*cur_Node)).Child_Right;
		}
	}
}
void Increas_No_Visit( Cluster_Node *cur_Node)
{
	Cluster_Node *tmp_Node = cur_Node;
	Max_No_Visit = (*tmp_Node).no_Visit + 1 > Max_No_Visit ?  (*tmp_Node).no_Visit + 1 : Max_No_Visit;
	while(tmp_Node != NULL)
	{
		(*tmp_Node).no_Visit++;
		tmp_Node = (*tmp_Node).Parent;
	}
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
//void NonRevisitingScheme_UnConstraint_Update(NonRevisitingScheme *cur_NRS, double *cur_Solution, double cur_Goodness)

int RecordRevisiting(double *cur_Chromosome, double cur_Chromonsome_Fitness,
		int no_gene, NonRevisitingScheme *cur_NRS, int *Current_ID, double **cur_Interval, vector<Selection_Strategy> &Select_Vector)
{
	double node_difference;
	Cluster_Node *cur_Node;
	int visitTimes = 0;
	RecordRevisiting_NodeSearch(cur_Chromosome, no_gene, &(*cur_NRS).Solution_Tree, &cur_Node);
	node_difference = Cluster_Node_Difference(cur_Node, cur_Chromosome, no_gene);
	if(node_difference == 0)
	{// the current node exist
	//	(*cur_Node).no_Visit++;
		Increas_No_Visit(cur_Node);
		visitTimes = (*cur_Node).no_Visit;

		// update, and store the best fitness value
		if (cur_Chromonsome_Fitness > (*cur_Node).Optimal_Fitness)
		{
			(*cur_Node).Optimal_Fitness = cur_Chromonsome_Fitness;
		}
	}else
	{//	insert data
		Cluster_Node *leaf_Node = cur_Node;
		double *cur_Data = cur_Chromosome;
		int i;
		int cur_Dimension = no_gene;
		double cur_Fitness = cur_Chromonsome_Fitness;
		double interval_diff;

		//cur_Node is the parent node that need to insert 
		(*cur_NRS).Solution_Tree.no_Leaf++;
		if(cur_Node == NULL)
			printf("NULL\n");
		// current node is root, then insert and replace it
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
		//note root node
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
			NonRevisitingScheme_UnConstraint_TreePath(leaf_Node, (*cur_NRS).Solution_Tree.cur_Interval, &(*cur_NRS).Solution_Tree.no_Leaf, (*cur_NRS).Solution_Tree.no_Dimension);
			visitTimes = 1;
			Archive_Leaf(leaf_Node, Select_Vector);
		}//	if((*leaf_Node).no_Visit == 0 && (*leaf_Node).Parent == NULL)
		if((*cur_NRS).Max_no_Leaf < (*cur_NRS).Solution_Tree.no_Leaf)
			(*cur_NRS).Max_no_Leaf = (*cur_NRS).Solution_Tree.no_Leaf;
	}//	if(node_difference == 0)
//	return node_difference;
	return visitTimes;
}

void Record_Cluster_Node_Write(Cluster_Node *cur_Node, double ***cur_Interval, int no_dimension, FILE *tree_ptr)
{
	int i;
	double min_bound, max_bound;
	if((*cur_Node).Child_Left != NULL || (*cur_Node).Child_Right != NULL)
	{
		min_bound = (*cur_Interval)[(*cur_Node).pseudo_Dimension][0];
		max_bound = (*cur_Interval)[(*cur_Node).pseudo_Dimension][1];

		if((*cur_Node).Child_Left != NULL)
		{
			(*cur_Interval)[(*cur_Node).pseudo_Dimension][0] = min_bound;
			(*cur_Interval)[(*cur_Node).pseudo_Dimension][1] = (*cur_Node).Threshold;
			Cluster_Node_Write((*cur_Node).Child_Left, cur_Interval, no_dimension, tree_ptr);
		}
		if((*cur_Node).Child_Right != NULL)
		{
			(*cur_Interval)[(*cur_Node).pseudo_Dimension][0] = (*cur_Node).Threshold;
			(*cur_Interval)[(*cur_Node).pseudo_Dimension][1] = max_bound;
			Cluster_Node_Write((*cur_Node).Child_Right, cur_Interval, no_dimension, tree_ptr);
		}
		(*cur_Interval)[(*cur_Node).pseudo_Dimension][0] = min_bound;
		(*cur_Interval)[(*cur_Node).pseudo_Dimension][1] = max_bound;

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
		fprintf(tree_ptr, "%1.0lf %d %d (%d. %d)\n", (*cur_Node).no_Visit, (*cur_Node).ID_Self, (*cur_Node).ID_Parent, (*cur_Node).Dimension, (*cur_Node).pseudo_Dimension);

	}
	return;
}

//Function to print out the Tree
void Record_Cluster_Tree_Write(Cluster_Tree *cur_Tree, char *tree_filename)
{
	FILE *tree_ptr;

	tree_ptr = fopen(tree_filename, "w");

	fprintf(tree_ptr, "%d\n", (*cur_Tree).no_Dimension);

	Cluster_Tree_SearchInterval_Set(cur_Tree);

	Record_Cluster_Node_Write((*cur_Tree).Root, &(*cur_Tree).cur_Interval, (*cur_Tree).no_Dimension, tree_ptr);

	fclose(tree_ptr);

	return;
}
void preOrder_Write_Into_Fstream(Cluster_Node *cur_Node, double ***cur_Interval,
		ostream &file_out, ostream &file_out_leaf, int dim, int axis_range, double *LBound, double *UBound)
{
	double min_bound, max_bound;
	    if((*cur_Node).pseudo_Dimension != -1)
	    {
			min_bound = (*cur_Interval)[(*cur_Node).pseudo_Dimension][0];
			max_bound = (*cur_Interval)[(*cur_Node).pseudo_Dimension][1];
			file_out<<"["<<(*cur_Interval)[(*cur_Node).pseudo_Dimension][0]<<",\t"<<(*cur_Interval)[(*cur_Node).pseudo_Dimension][1]<<"]\t"
					<< setw(4) << setfill(' ')<<(*cur_Node).ID_Self<<">"
					<<setw(4) << setfill(' ')<<(*cur_Node).ID_Parent
					<<setw(10) << setfill(' ')<<(*cur_Node).no_Visit
					<<"\t("<<setw(2)<<  setfill(' ')<<(*cur_Node).pseudo_Dimension<<",\t"
					<<setw(2)<<  setfill(' ')<<(*cur_Node).Dimension<<")\t";
			file_out<<(*cur_Node).Optimal_Fitness<<"\t:\t";
			for(int j= 0; j < dim; j++)
				file_out<<(*cur_Node).Optimal_Data[j]<<"\t";
			file_out<<"|\t";
			for(int j= 0; j < dim; j++)
			{
				double tmp_pop = (*cur_Node).Optimal_Data[j] / (double)(axis_range - 1) * (UBound[j] - LBound[j]) + LBound[j];
				file_out<<tmp_pop<<"\t";
			}
			file_out<<endl;
			if((*cur_Node).Child_Left != NULL)
			{
				(*cur_Interval)[(*cur_Node).pseudo_Dimension][0] = min_bound;
				(*cur_Interval)[(*cur_Node).pseudo_Dimension][1] = (*cur_Node).Threshold;
				preOrder_Write_Into_Fstream((*cur_Node).Child_Left, cur_Interval, file_out, file_out_leaf, dim, axis_range, LBound, UBound);
			}

			if((*cur_Node).Child_Right!=NULL)
			{
				(*cur_Interval)[(*cur_Node).pseudo_Dimension][0] = (*cur_Node).Threshold;
				(*cur_Interval)[(*cur_Node).pseudo_Dimension][1] = max_bound;
				preOrder_Write_Into_Fstream((*cur_Node).Child_Right, cur_Interval, file_out, file_out_leaf, dim, axis_range, LBound, UBound);
			}
			(*cur_Interval)[(*cur_Node).pseudo_Dimension][0] = min_bound;
			(*cur_Interval)[(*cur_Node).pseudo_Dimension][1] = max_bound;
		}else
			if((*cur_Node).pseudo_Dimension == -1)
			{
				for(int j = 0; j < dim; j++)
				{
					file_out_leaf<<"["<<(*cur_Interval)[j][0]<<",\t"<<(*cur_Interval)[j][1]<<"]\t";
				}
				file_out_leaf<< setw(4) << setfill(' ')<<(*cur_Node).ID_Self<<"\t>\t"
						<<setw(4) << setfill(' ')<<(*cur_Node).ID_Parent
						<<setw(10) << setfill(' ')<<(*cur_Node).no_Visit
						<<"\t("<<setw(2)<<  setfill(' ')<<(*cur_Node).pseudo_Dimension<<",\t"
						<<setw(2)<<  setfill(' ')<<(*cur_Node).Dimension<<")\t";
				file_out_leaf<<(*cur_Node).Optimal_Fitness<<"\t:\t";
				for(int j= 0; j < dim; j++)
					file_out_leaf<<(*cur_Node).Optimal_Data[j]<<"\t";
				file_out_leaf<<"|\t";
				for(int j= 0; j < dim; j++)
				{
					double tmp_pop = (*cur_Node).Optimal_Data[j] / (double)(axis_range - 1) * (UBound[j] - LBound[j]) + LBound[j];
					file_out_leaf<<tmp_pop<<"\t";
				}
				file_out_leaf<<endl;
			}
}
void preOrder_Write_Into_File(Cluster_Tree *cur_Tree, double ***cur_Interval, const char *filename,
		const char *filename_leaf,	int func, int run, int maxLeaf, int dim, int axis_range, double *LBound, double *UBound)
{
	std::ofstream File_write_out(filename);
	std::ofstream File_write_out_leaf(filename_leaf);
	File_write_out<<maxLeaf<<endl;
	File_write_out_leaf<<"Leaf Node:\t"<<maxLeaf<<endl;
	File_write_out_leaf<<"Node UsedAfterIni:\t"<<No_Leaf_Initialization<<endl;
	File_write_out_leaf<<"Mean Probability:\t"<<Mutate_Probability<<endl;
	File_write_out_leaf<<"Mean AxisRange:\t"<<Mutate_On_Axis_Range<<endl;
	File_write_out_leaf<<"Mean NearestSpace:\t"<<Mutate_To_Nearest_Space<<endl;
	File_write_out_leaf<<"Max NoVisit:\t"<<Max_No_Visit<<endl;
	
	/*File_write_out_leaf << maxLeaf << endl;
	File_write_out_leaf << No_Leaf_Initialization << endl;
	File_write_out_leaf << Mutate_Probability << endl;
	File_write_out_leaf << Mutate_On_Axis_Range << endl;
	File_write_out_leaf <<  Mutate_To_Nearest_Space << endl;
	File_write_out_leaf << Max_No_Visit << endl;*/

	Cluster_Tree_SearchInterval_Set(cur_Tree);
	preOrder_Write_Into_Fstream((*cur_Tree).Root, cur_Interval, File_write_out, File_write_out_leaf, dim, axis_range, LBound, UBound);
	File_write_out_leaf.close();
	File_write_out.close();
}

void preOrder_Write_Into_Fstream(Cluster_Node *cur_Node, double ***cur_Interval,
	ostream &file_out_leaf, int dim, int axis_range, double *LBound, double *UBound)
{
	double min_bound, max_bound;
	if ((*cur_Node).pseudo_Dimension != -1)
	{
		min_bound = (*cur_Interval)[(*cur_Node).pseudo_Dimension][0];
		max_bound = (*cur_Interval)[(*cur_Node).pseudo_Dimension][1];
		if ((*cur_Node).Child_Left != NULL)
		{
			(*cur_Interval)[(*cur_Node).pseudo_Dimension][0] = min_bound;
			(*cur_Interval)[(*cur_Node).pseudo_Dimension][1] = (*cur_Node).Threshold;
			preOrder_Write_Into_Fstream((*cur_Node).Child_Left, cur_Interval, file_out_leaf, dim, axis_range, LBound, UBound);
		}

		if ((*cur_Node).Child_Right != NULL)
		{
			(*cur_Interval)[(*cur_Node).pseudo_Dimension][0] = (*cur_Node).Threshold;
			(*cur_Interval)[(*cur_Node).pseudo_Dimension][1] = max_bound;
			preOrder_Write_Into_Fstream((*cur_Node).Child_Right, cur_Interval, file_out_leaf, dim, axis_range, LBound, UBound);
		}
		(*cur_Interval)[(*cur_Node).pseudo_Dimension][0] = min_bound;
		(*cur_Interval)[(*cur_Node).pseudo_Dimension][1] = max_bound;
	}
	else
		if ((*cur_Node).pseudo_Dimension == -1)
		{
			for (int j = 0; j < dim; j++)
			{
				file_out_leaf << "[" << (*cur_Interval)[j][0] << ",\t" << (*cur_Interval)[j][1] << "]\t";
			}
			file_out_leaf << setw(4) << setfill(' ') << (*cur_Node).ID_Self << "\t>\t"
				<< setw(4) << setfill(' ') << (*cur_Node).ID_Parent
				<< setw(10) << setfill(' ') << (*cur_Node).no_Visit
				<< "\t(" << setw(2) << setfill(' ') << (*cur_Node).pseudo_Dimension << ",\t"
				<< setw(2) << setfill(' ') << (*cur_Node).Dimension << ")\t";
			file_out_leaf << (*cur_Node).Optimal_Fitness << "\t:\t";
			for (int j = 0; j < dim; j++)
				file_out_leaf << (*cur_Node).Optimal_Data[j] << "\t";
			file_out_leaf << "|\t";
			for (int j = 0; j < dim; j++)
			{
				double tmp_pop = (*cur_Node).Optimal_Data[j] / (double)(axis_range - 1) * (UBound[j] - LBound[j]) + LBound[j];
				file_out_leaf << tmp_pop << "\t";
			}
			file_out_leaf << endl;
		}
}
void preOrder_Write_Into_File(Cluster_Tree *cur_Tree, double ***cur_Interval, 
	const char *filename_leaf, int func, int run, int maxLeaf, int dim, int axis_range, double *LBound, double *UBound)
{
	std::ofstream File_write_out_leaf(filename_leaf);
	File_write_out_leaf << "Leaf Node:\t" << maxLeaf << endl;
	File_write_out_leaf << "Node UsedAfterIni:\t" << No_Leaf_Initialization << endl;
	File_write_out_leaf << "Mean Probability:\t" << Mutate_Probability << endl;
	File_write_out_leaf << "Mean AxisRange:\t" << Mutate_On_Axis_Range << endl;
	File_write_out_leaf << "Mean NearestSpace:\t" << Mutate_To_Nearest_Space << endl;
	File_write_out_leaf << "Max NoVisit:\t" << Max_No_Visit << endl;

	/*File_write_out_leaf << maxLeaf << endl;
	File_write_out_leaf << No_Leaf_Initialization << endl;
	File_write_out_leaf << Mutate_Probability << endl;
	File_write_out_leaf << Mutate_On_Axis_Range << endl;
	File_write_out_leaf <<  Mutate_To_Nearest_Space << endl;
	File_write_out_leaf << Max_No_Visit << endl;*/

	Cluster_Tree_SearchInterval_Set(cur_Tree);
	preOrder_Write_Into_Fstream((*cur_Tree).Root, cur_Interval, File_write_out_leaf, dim, axis_range, LBound, UBound);
	File_write_out_leaf.close();
}