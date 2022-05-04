#ifndef OPTIMIZATION_TESTFUNCTION_KERNEL

#define OPTIMIZATION_TESTFUNCTION_KERNEL

// global variables
extern int Optimization_TestFunction_Dimension;
extern double Optimization_TestFunction_AxisRange;

extern int Optimization_TestFunction_RealRoyalRoadTau;
extern int Optimization_TestFunction_MultiRealValue_BlockSize;

extern double Optimization_TestFunction_LinearWeigth[256];
extern double Optimization_TestFunction_XORTReeNode[2][256];
extern double *Optimization_TestFunction_TruthTableValue;

extern double Optimization_TestFunction_Unbounded_Bias;

extern double Optimization_TestFunction_NoiseLevel;
extern double (*Optimization_TestFunction_TruthTableFunction)(double*);
extern double (*Optimization_TestFunction_Pointer)(double*);


// variables
extern int Optimization_TestFunction_RotationFlag;
extern double **Optimization_TestFunction_RotationMatrix, *Scaled_Solution, *Rotated_Solution;

void   Optimization_TestFunction_RotationMatrix_Construction();
void   Optimization_TestFunction_RotationMatrix_Load();
void   Optimization_TestFunction_RotationMatrix_Destruction();
void   Optimization_TestFunction_RotationSolution(double*, double, double);

void	Optimization_TestFunction_Manager_RealValued(int, char*);
void	Optimization_TestFunction_Manager_PesudoBoolean(int, char*);
void	Optimization_TestFunction_Information(double**, double**, double**, double*, int);


//========= Uni-modal Function
double Optimization_TestFunction_Planar(double*);
double Optimization_TestFunction_Sphereical(double*);
double Optimization_TestFunction_Schwefel2_22(double*);
double Optimization_TestFunction_Schwefel1_02(double*);
double Optimization_TestFunction_Schwefel2_21(double*);
double Optimization_TestFunction_Rosenbrock(double*);
double Optimization_TestFunction_Quartic(double*);

//========= Multi-modal Function
double Optimization_TestFunction_GeneralizedRastrigin(double*);
double Optimization_TestFunction_GeneralizedGriewank(double*);
double Optimization_TestFunction_Schwefel2_26(double*);
double Optimization_TestFunction_Ackley(double*);
double Optimization_TestFunction_Foxholes(double*);
double Optimization_TestFunction_SixHumpCamelBack(double*);
double Optimization_TestFunction_Branin(double*);
double Optimization_TestFunction_GoldsteinPrice(double*);
double Optimization_TestFunction_HighConditionedElliptic(double*);
double Optimization_TestFunction_Weierstrass(double*);
double Optimization_TestFunction_Sinc(double*);
double Optimization_TestFunction_Levy(double*);
double Optimization_TestFunction_Pern(double*);
double Optimization_TestFunction_Zakharov(double*);
double Optimization_TestFunction_Alpine(double*);
double Optimization_TestFunction_Pathological(double*);
double Optimization_TestFunction_InvertedCosineWave(double*);
double Optimization_TestFunction_InvertedCosineMixture(double*);
double Optimization_TestFunction_EpistaticMichalewicz(double*);
double Optimization_TestFunction_LevyMontalvo(double*);
double Optimization_TestFunction_Neumaier3(double*);
double Optimization_TestFunction_OddSquare(double*);
double Optimization_TestFunction_Paviani(double*);
double Optimization_TestFunction_Periodic(double*);
double Optimization_TestFunction_Salomon(double*);
double Optimization_TestFunction_Shubert(double*);
double Optimization_TestFunction_Sinusoidal(double*);
double Optimization_TestFunction_Michalewicz(double*);
double Optimization_TestFunction_Whitely(double*);
double Optimization_TestFunction_Whitely(double*);
double Optimization_TestFunction_Easom(double*);

//========= Hybrid Composition Function
double Optimization_TestFunction_HybridComposition01(double*);

//========= Pesudo Boolean Function
double Optimization_TestFunction_MaxOne(double*);
double Optimization_TestFunction_RealRoyalRoad(double*);
double Optimization_TestFunction_Linear(double*);
double Optimization_TestFunction_AlmostPositive(double*);
double Optimization_TestFunction_InversePseudoModular(double*);
double Optimization_TestFunction_RealValue(double*);
double Optimization_TestFunction_ShortPathConstant(double*);
double Optimization_TestFunction_XORTree(double*);
double Optimization_TestFunction_MultiRealValue(double*);
double Optimization_TestFunction_TruthTable(double*);
double Optimization_TestFunction_InverseRealRoyalRoad2001(double*);

void   Optimization_TestFunction_TruthTableLoad(char*);
void   Optimization_TestFunction_TruthTableConvert(double (*fn_ptr)(double*), int);
void   Optimization_TestFunction_TruthTableRandomSwapping();
void   Optimization_TestFunction_TruthTableDestruction();

#endif
