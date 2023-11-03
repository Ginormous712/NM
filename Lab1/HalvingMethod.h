#pragma once
#include <iostream>
#include <cmath>

using namespace std;
// Given function  x^3 - 4x^2 - 15x + 18 = 0
double f(double x)
{
    return x * x * x - 4 * x * x - 15 * x + 18;
}

double halvingMethod(double a, double b, double (*f)(double), double epsilon)
{
    double c; 
    int n = 1;
    
    cout << endl << "**************************" << endl;
    cout << "Halving Method" << endl;
    cout << "**************************" << endl;
    if (f(a) * f(b) > 0) { 
        std::cout<<("Error! Can't use method\n"); 
        return 0; 
    }

    std::cout << "Number od iterations:\t";
    std::cout << (int)(log((b - a) / epsilon) / log(2) + 1) << std::endl;

    while (abs(b - a) > epsilon) {
        c = (a + b) / 2.0;
        double fc = f(c);
        double fa = f(a);
        double fb = f(b);

        std::cout << n << "|\ta:\t" << a << "\tc:\t" << c << "\tb:\t" << b << std::endl;
        if (fa * fc > 0) { a = c;}
        else { b = c;}
        n++;
    }
    return c;
}
