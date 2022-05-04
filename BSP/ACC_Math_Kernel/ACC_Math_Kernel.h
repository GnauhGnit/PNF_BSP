#ifndef ACC_MATH_KERNEL

#define ACC_MATH_KERNEL

#define INF 99999999999.0	
#define PI 3.14159265

#define EXPONENTIAL_NO_DIVID 5096
#define NORMALIZE_SQUARE_NO_DIVID 5096

//double Exponential_Width_neg, Exponential_Width_pos, *Exponential_Array;
//double Normalize_Square_Width, *Normalize_Square_Array;

extern double Exponential_Width_neg, Exponential_Width_pos, *Exponential_Array;
extern double Normalize_Square_Width, *Normalize_Square_Array;

void Exponential_Construction();
void Exponential_Destruction();
double Exponential(double);

void Normalize_Square_Construction();
void Normalize_Square_Destruction();
double Normalize_Square(double);

double randu();
double randn(double, double);

void Sorting_Indexed(double**, int**, int);
void Sorting_UnIndexed(double**, int);


double IMGau_Gau(double*, double*, double, double, int, double**);
double IMGau(double*, double, int, double**);
double IGau_Gau(double, double, double, double, double, double);
double IGau(double, double, double, double);

int	   Matrix_Inversion(double**, double**, int, double*);
void   OrthogonalMatrix_Generator(int, double***);
void   Orthonormal_Basis(double*, int, double***);

double Min(double*, int);
double Max(double*, int);

double Sigma(double, double);
double dSigma(double, double);

void Dec2Bin(int, double**, int);
double nChoosek(int, int);

int Sign(double);

void Pause();



#endif
