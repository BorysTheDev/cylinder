#pragma once
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

	
	double GetConstEps() const;
    int32_t GetSize() const;
    int32_t GetTestPointsNum() const;
	
    friend std::istream& operator>>(std::istream& input, Given& given);

private:
    int32_t N;
    int32_t TestPointsNum;
};

double Given::GetConstEps() const
{
    return M_PI / eps;
}

int32_t Given::GetSize() const
{
    return N - 2;
}

int32_t Given::GetTestPointsNum() const
{
    return TestPointsNum;
}

std::istream& operator>>(std::istream& input, Given& given)
{
    double real, img;
    input >> given.eps >> real >> img >> given.N >> given.TestPointsNum;
    given.kappa = cmpx(real, img);
    return input;
}

