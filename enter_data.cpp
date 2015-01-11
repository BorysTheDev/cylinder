#include "enter_data.h"
#include "chebyshev.h"
#include "hankel.h"
#include "simple_math.h"

using namespace mth;
void Enter_Data ( cmpx kappa, double ConstEps, Matrix& M)
{
    // N: Integer:= M.all'Last(1);
    const cmpx A0 = 1.0 / (2.0*kappa);
    const cmpx ConstA1 = -kappa / 4.0;
    const cmpx ConstG0 = ((M_PI * kappa) / cmpx(0, 2.0)) * j0(kappa.real()) * h1(kappa.real());

    double SumCos = 0.0;

    for (int i = 1; i <= 2000; i++)  // what is 2000?
    {
        SumCos += 2.0 / sqr(M_PI * i);
    }

    int32_t N = M.size() + 1;
    double dn = N;

    for (int tempI = 1; tempI < N; tempI++)
    {
        for (int tempJ = 1; tempJ < N; tempJ++)
        {
            if (tempI == tempJ)
            {

                M[tempI][tempJ] = -(A0 / (sqr(M_PI - ConstEps))) * (dn / 2.0) +
                        ( (1.0 - sqr(uchebNodes<double>(N, tempI)) ) / dn ) *
                        (-ConstA1*(log(2.0)+ 2.0 * Sum_Base_Elem(tempI,tempI,N) + (1.0 / (2.0 * dn)) )
                         + (A0 / 4.0) * SumCos + (ConstA1 * log(abs(M_PI - ConstEps)) - ConstG0/2.0)-
                         KernelK(kappa,A0,ConstA1,ConstEps,tempI,tempI,N) );

            }
            else
            {
                double sign = (tempI + tempJ) % 2 ? -1.0 : 1;
                M[tempI][tempJ] = (
                            ( 1.0 - sqr(uchebNodes<double>(tempJ, N)) ) / dn *
                            ( 1.0 - (sign)/
                               sqr(uchebNodes<double>(tempI, N) - uchebNodes<double>(tempJ,N)))*
                              (A0 / sqr(M_PI - ConstEps) ) -
                              ConstA1*(log(2.0) + 2.0 * Sum_Base_Elem(tempI, tempJ, N) +
                                       sign / double(2*N) ) +
                              ( A0 / 2.0 )*(1.0 / (2.0 * sqr(sin(((M_PI - ConstEps) / 2.0)*
                                                      (uchebNodes<double>(tempI,N) - uchebNodes<double>(tempJ,N))))) -
                                        2.0 / (sqr(M_PI - ConstEps) * sqr(uchebNodes<double>(tempI,N) -
                                                                uchebNodes<double>(tempJ,N))))+
                              ConstA1 * log(abs((sin(((M_PI-ConstEps)/2.0)*
                                                   (uchebNodes<double>(tempI,N)-uchebNodes<double>(tempJ,N))))/
                                              (((M_PI-ConstEps)/2.0)*(uchebNodes<double>(tempI,N)-uchebNodes<double>(tempJ,N)))))+
                              ConstA1*log(abs(M_PI - ConstEps))-ConstG0/2.0-
                              KernelK(kappa, A0, ConstA1, ConstEps, tempJ, tempI, N));
            }


        }
    }

}

cmpx Fun_Gn(cmpx kappa, int Nu)
{
    double dn = Nu;
    return ( (M_PI*kappa) / cmpx(0,2) ) * ( -_jn(Nu+1, kappa.real() )+ (dn/kappa) *
                                            _jn(Nu, kappa.real())) * (-h1(Nu+1, kappa.real()) + (dn/kappa) * h1(Nu, kappa.real()));


}

cmpx KernelK(cmpx kappa, cmpx ConstA0, cmpx ConstA1, double ConstEps, int K, int J, int N)
{
    const int tempN = 150;
    cmpx Tempsum = 0;
    for (int tempI = 1; tempI <= tempN; tempI++ )
    {
        Tempsum += ( Fun_Gn(kappa, tempI) - ConstA0 * double(tempI)
                     - ConstA1 * (1.0 / tempI)) * cos(tempI * (M_PI - ConstEps)
                                                      * (uchebNodes<double>(J, N) - uchebNodes<double>(K, N)));
    }

    return Tempsum;
}

double Pol_Cheb_I(int K, int tempN, int N)
{
    return cos( tempN * acos(uchebNodes<double>(K,N)) );
}

double Sum_Base_Elem (int K, int J, int N)
{
    double Tempsum =0.0;
    for (int Tempr = 1; Tempr < N - 1; Tempr++)
    {
        Tempsum += (Pol_Cheb_I(K,Tempr,N) * Pol_Cheb_I(J,Tempr,N)) / Tempr;
    }
    return Tempsum;
}
