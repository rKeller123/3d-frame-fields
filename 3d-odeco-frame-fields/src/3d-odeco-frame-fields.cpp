// 3d-odeco-frame-fields.cpp : Defines the entry point for the application.

#include "3d-odeco-frame-fields.h"

using namespace std;

int main()
{
	double l_x = 1.0;
	double l_y = 1.0;
	double l_z = 1.0;

	double alpha = 1.0;
	double beta  = 2.0;
	double gamma = 3.0;

	vec15 f = odeco_frame(alpha, beta, gamma, l_x, l_y, l_z);

	vec15 y = vec15::Constant(1);

	cout << setprecision(17);

	cout << "Odeco frame: " << endl << f << endl << endl;

	cout << "Objective function: " << objective_function(y, alpha, beta, gamma, l_x, l_y, l_z) << endl << endl;

	cout << "Gradient: " << endl << gradient(y, alpha, beta, gamma, l_x, l_y, l_z) << endl << endl;

	mat6 H = hessian(y, alpha, beta, gamma, l_x, l_y, l_z);

	cout << "Hessian: " << endl << H << endl << endl;

	cout << "Hessian Eigenvalues: " << endl << H.eigenvalues() << endl << endl;

	cout << setprecision(5);
	cout << "Hessian - Hessian^transpose: " << endl << H - H.transpose() << endl << endl;


	return 0;
}
