/*
 * RecordRevisitingScheme_Kernel.h
 *
 *  Created on: Jul 29, 2017
 *      Author: hting
 */

#ifndef BSP_RECORDREVISITINGSCHEME_KERNEL_RECORDREVISITINGSCHEME_KERNEL_H_
#define BSP_RECORDREVISITINGSCHEME_KERNEL_RECORDREVISITINGSCHEME_KERNEL_H_

#include "../../BSP/ACC_Math_Kernel/ACC_Math_Kernel.h"
#include "../../BSP/ClusterTree_Kernel/ClusterTree_Kernel.h"
#include "../../BSP/GeneticAlgorithm_Kernel/GeneticAlgorithm_Kernel.h"
#include "../../BSP/NonRevisitingScheme_Kernel/NonRevisitingScheme_Kernel.h"
#include "../../BSP/SelectionStrategey_Kernel/Selection_Strategy.h"
#include <fstream>
#include <vector>

//int RecordRevisiting(double *cur_Chromosome, double cur_Chromonsome_Fitness,
//		int no_gene, NonRevisitingScheme *cur_NRS, int *Current_ID, double **cur_Interval);
int RecordRevisiting(double *cur_Chromosome, double cur_Chromonsome_Fitness,
	int no_gene, NonRevisitingScheme *cur_NRS, int *Current_ID, double **cur_Interval, std::vector<Selection_Strategy> &Select_Vector);


void Record_Cluster_Node_Write(Cluster_Node *cur_Node, double ***cur_Interval, int no_dimension, FILE *tree_ptr);

void Record_Cluster_Tree_Write(Cluster_Tree *cur_Tree, char *tree_filename);
//void preOrder_2(Cluster_Node *cur_Node);
//void preOrder(Cluster_Node *cur_Node, double ***cur_Interval);
//
//void preOrder_2(Cluster_Node *cur_Node, char *tree_filename);
//void preOrder(Cluster_Node *cur_Node, double ***cur_Interval, char *tree_filename);

//new file
void preOrder_Write_Into_File(Cluster_Tree *cur_Tree, double ***cur_Interval, const char *filename,
		const char *filename_leaf,	int func, int run, int maxLeaf, int dim, int axis_range, double *LBound, double *UBound);
void preOrder_Write_Into_File(Cluster_Tree *cur_Tree, double ***cur_Interval,
	const char *filename_leaf, int func, int run, int maxLeaf, int dim, int axis_range, double *LBound, double *UBound);
#endif /* BSP_RECORDREVISITINGSCHEME_KERNEL_RECORDREVISITINGSCHEME_KERNEL_H_ */
