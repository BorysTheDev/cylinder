#pragma once
#include <vector>
#include "given.h"


typedef std::vector<std::vector<cmpx>> Matrix;

cmpx KernelK(cmpx kappa, cmpx ConstA0, cmpx ConstA1, double ConstEps, int K, int J, int N);

void Enter_Data ( cmpx kappa, double ConstEps, Matrix& M);

cmpx Fun_Gn(cmpx kappa, int Nu);

cmpx KernelK(cmpx kappa, cmpx ConstA0, cmpx ConstA1, double ConstEps, int K, int J, int N);

double Pol_Cheb_I(int K, int TempN, int N);

double Sum_Base_Elem (int K, int J, int N);

//Fun_T_0 is mth::uchebNodes<double>
