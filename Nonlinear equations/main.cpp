#include <iostream>
#include "HalvingMethod.h"
//#include "NewtonMethod.h"
//#include "FixedPointIteration.h"

using namespace std;

double epsilon = 0.001;

int main() {

	//cout << "Enter epsilon (0.001 by default):\n";
	//cin >> epsilon;
	cout << "Root is: " << halvingMethod(1, 6, f, epsilon);
	//cout << "Root is: " << newtonMethod(3.5, 5.5, 4.5, f, epsilon);
	//cout << "Root is: " << fixedPointMethod(0.5, 1.5, 0.5, f, g, epsilon);
	return 0;
}