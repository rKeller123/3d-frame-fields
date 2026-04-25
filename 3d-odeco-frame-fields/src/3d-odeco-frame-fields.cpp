// 3d-odeco-frame-fields.cpp : Defines the entry point for the application.

#include "3d-odeco-frame-fields.h"

using namespace std;

int main()
{
	double l_x = 1.0;
	double l_y = 1.0;
	double l_z = 1.0;

	double alpha = 0.0;
	double beta  = 0.0;
	double gamma = 0.0;

	vec15 f = odeco_frame(alpha, beta, gamma, l_x, l_y, l_z);

	cout << setprecision(17);

	cout << objective_function(vec15::Constant(1), alpha, beta, gamma, l_x, l_y, l_z) << endl;

	return 0;
}
