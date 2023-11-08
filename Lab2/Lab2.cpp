#include <iostream>
#include "Matrix.h"

#define SPACE "*************************************************\n"

using namespace std;



int main()
{
    vector<vector<double>> data = {
        {4.0, 3.0, 1.0, 0.0},
        {-2.0, 2.0, 6.0, 1.0},
        {0.0, 5.0, 2.0, 3.0},
        {0.0, 1.0, 2.0, 7.0}
    };

    vector<double> b = { 14,31,33,45 };

    try 
    {
        SquareMatrix a(data);
        cout << SPACE;
        cout << "The square matrix was created successfully." << endl;
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
        cout << "Vector b:\t";
        
        for (size_t i = 0; i < b.size(); i++)
        {
            cout << b[i] << " ";
        }
        cout << endl;

        vector<double> result = a.solveLinearSystemGauss(b);
        cout << "Solution:\t";
        for (size_t i = 0; i < b.size(); i++)
        {
            cout << "x" << i + 1  << " = " << result[i] << ";\t";
        }
        cout << endl;
    }
    catch (const invalid_argument& e) 
    {
        cerr << "Error: " << e.what() << endl;
    }

}