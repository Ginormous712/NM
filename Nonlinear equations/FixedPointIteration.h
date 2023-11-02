#pragma once
#include <iostream>
#include <cmath>
using namespace std;

// Given function  x^3 - 7x^2 - x + 7 = 0
double f(double x)
{
    return x * x * x - 7 * x * x - x + 7;
}

double g(double x)
{
    return sqrt((x*x*x - x + 7)/7);
}

// Function to find the root
double fixedPointMethod(double a, double b, double x0, 
    double (*f)(double), double (*g)(double), double epsilon)
{
    cout << "**************************" << endl;
    cout << "Fixed-Point Method" << endl;
    cout << "**************************" << endl;
    
    if (f(a) * f(b) > 0) {
        std::cout << ("Error! Can't use method\n");
        return 0;
    }
    double q = 0.365;
    std::cout << "Number od iterations:\t";
    std::cout << (int)(log((abs(g(x0) - x0)) / ((1 - q) * epsilon)) / log(1/q)) + 1 << std::endl;


    double x;
    int n = 1;
    do
    {
        x = g(x0);
        std::cout << n << "|\tx:\t" << x << "\t x - g(x)\t" << x - g(x) << std::endl;
        x0 = x;
        n++;
        
    } while (abs((x - g(x)) > epsilon));

    //std::cout << "The value of the root is : " << x;
    return x;
}