
#define Message  printf
#define Error(s) {printf(s);exit(0);}

#define NUM 100  //The number of AD;
#define PI 3.1415926
#define eps 1E-15

//double X_PMT[NUM];
//double Y_PMT[NUM];
//double Z_PMT[NUM];  //The position of PMT;
//double Q_expected[NUM];  //The expected charge; 

//double Grad[4];  //The gradient;

double Likelihood_v;  //The likelihood value;

const int ReadDataFromFile = 1;  //If read data from file or random real estate;
const double t = 1E-4;  //The step length;

//Other parameters/////
//double W;  //Other parameter's ji;
const double Lambda_a = 11.0;
const double Effi_QE = 1.0;
const double Light_yield = 6.489E-3;
const double S_pmt = 0.0025;
