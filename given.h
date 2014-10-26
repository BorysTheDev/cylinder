#include <complex>
#include <cstdint>
#include <istream>
#include <cmath>

typedef std::complex<double> cmpx;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

struct Given
{
	double eps;
    cmpx kappa;
	int32_t N;
	
	double GetConstEps() const;
	
    friend std::istream& operator>>(std::istream& input, Given& given);
};

double Given::GetConstEps() const
{
    return M_PI / eps;
}

std::istream& operator>>(std::istream& input, Given& given)
{
    double real, img;
    input >> given.eps >> real >> img >> given.N;
    given.kappa = cmpx(real, img);
    return input;
}

