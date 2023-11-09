#include <iostream>
#include "Matrix.h"

#define SPACE "*************************************************\n"

using namespace std;

int main()
{
    vector<vector<double>> data1 = {
        {4.0, 3.0, 1.0, 0.0},
        {-2.0, 2.0, 6.0, 1.0},
        {0.0, 5.0, 2.0, 3.0},
        {0.0, 1.0, 2.0, 7.0}
    };

    vector<vector<double>> data2 = {
        {1.0, 2.0, 0.0},
        {2.0, 2.0, 4.0},
        {0.0, 4.0, 3.0}
    };

    vector<vector<double>> data3 = {
        {4.0, 0.0, 1.0, 1.0},
        {0.0, 3.0, 0.0, 1.0},
        {1.0, 0.0, 2.0, 0.0},
        {1.0, 1.0, 0.0, 5.0}
    };

    vector<double> b1 = { 14,31,33,45 };
    vector<double> b2 = { 11,34,31 };
    vector<double> b3 = { 11,10,7,23 };

    try 
    {
        cout << SPACE;
        cout << "Gaussian method" << endl;
        SquareMatrix a(data1);
        cout << SPACE;

        cout << "Matrix A:\n" << a;
        cout << SPACE;

        cout << "Determiant of matrix A:\t" << a.determinant() << endl;
        cout << SPACE;

        cout << "Invesre matrix to A:\n";
        SquareMatrix inverse = a.findInverseMatrixGauss();
        cout << inverse;

        cout << SPACE;
        cout << "Solving linear system using Gaussian method:\n";
        vector<double> result1 = a.solveLinearSystemGauss(b1);
        a.solveAndPrint(b1, result1);

        cout << SPACE;
        cout << "Tridiagonal matrix algorithm" << endl;
        TridiagonalMatrix b(data2);
        cout << SPACE;
        cout << "Solving linear system using Tridiagonal matrix algorithm:\n";
        vector<double> result2 = b.tridiagonalSolve(b2);
        b.solveAndPrint(b2, result2);


        
    }
    catch (const invalid_argument& e) 
    {
        cerr << "Error: " << e.what() << endl;
    }

}