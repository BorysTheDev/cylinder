#include <iostream>
#include "given.h"
#include "gauss.h"
#include "enter_data.h"

int main()
{
    //---test examples------------------



    //--------------------------------------------------
    //unclear variables
    double K = 4.399;
    double A = 1;
    double empty = 0;
    std::cout << "Please, give us real Eps (>= 0), Kappa (real+i*real), positise integer N and test points number :\n";
    Given given;
	std::cin >> given;
	
    Matrix matrix(given.GetSize());
    for (auto& v : matrix)
    {
        v.resize(given.GetSize());
    }

    std::cout << " matrix was allocated\n";

    for (int i = 0; i < given.GetTestPointsNum(); i++)
    {
        cmpx my_kappa(K * A, empty);
        std::cout << "Det(" << K << ")\n";
    }

    return 0;
}
