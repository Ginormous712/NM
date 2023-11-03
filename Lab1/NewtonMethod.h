#pragma once
#include <iostream>
#include <cmath>

using namespace std;
// Given function  x^3 - 10x^2 + 11x + 70 = 0

bool max(double a, double b)
{
    if (a >= b)
    {
        return a;
    }
    else
    {
        return b;
    }
}

double f(double x)
{
    return x * x * x - 10 * x * x + 11 * x + 70;
}

double df(double x)
{
    return 3 * x * x - 20 * x + 11;
}

double ddf(double x)
{
    return 6 * x - 20;
}

// Function to find the root
double newtonMethod(double a, double b, double x0, double (*f)(double), double epsilon)
{
    cout << endl << "**************************" << endl;
    cout << "Newton's Method" << endl;
    cout << "**************************" << endl;
    if (f(a) * f(b) > 0) {
        std::cout << ("Error! Can't use method\n");
        return 0;
    }
    double m1 = abs(df(b));
    double m2 = abs(ddf(b));
    double delta = max(b-x0, x0-a);
    double q = m2 * delta / (2 * m1);
    std::cout << "Number od iterations:\t";
    std::cout << (int)(log(log(delta/epsilon)/log(1/q))/log(2) + 1) + 1 << std::endl;



    int n = 1;
    double x = x0;
    double h;
    do 
    {
        h = f(x) / df(x);
        std::cout << n << "|\tx:\t" << x << "\tf(x)\t" << f(x) << std::endl;
        x = x - h;
        n++;

    } while (abs(f(x)) > epsilon);
    std::cout << n << "|\tx:\t" << x << "\tf(x)\t" << f(x) << std::endl;
    //std::cout << "The value of the root is : " << x;
    return x;
}

