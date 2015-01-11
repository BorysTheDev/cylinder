#include "given.h"

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
