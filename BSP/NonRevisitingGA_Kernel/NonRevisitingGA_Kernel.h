#ifndef NONREVISITINGGA_KERNEL

#define NONREVISITINGGA_KERNEL

#ifndef CLUSTERTREE_KERNEL
	#include"../../BSP/ClusterTree_Kernel/ClusterTree_Kernel.h"
#endif

#ifndef GENETICALGORITHM_KERNEL
	#include"../../BSP/GeneticAlgorithm_Kernel/GeneticAlgorithm_Kernel.h"
#endif

#ifndef NONREVISITINGSCHEME_KERNEL
	#include"../../BSP/NonRevisitingScheme_Kernel/NonRevisitingScheme_Kernel.h"
#endif

typedef struct NONREVISITING_GA {

	GeneticAlgorithm GA_Info;		// The original GA

	NonRevisitingScheme NRS_Info;	// The non-revisiting Scheme

} NonRevisiting_GA;

#define NONREVISITINGGA_MAX_GENERATION 50000
#define NONREVISITINGGA_STOPPING_GENERATION 10

void NonRevisitingGA_Construction(NonRevisiting_GA*, int, int, int, double, double, int, int);
void NonRevisitingGA_Destruction(NonRevisiting_GA);

//============ UnConstraint Optimization

void NonRevisitingGA_UnConstraint_Standard(NonRevisiting_GA*, double, int, double, int, double (*fn_ptr)(double*));
void NonRevisitingGA_UnConstraint_Mutation(double**, int, NonRevisitingScheme*);
#endif
