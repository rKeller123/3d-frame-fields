#include "odeco.h"


using namespace std;

vec15 reference_frame(double l_x, double l_y, double l_z)
{
	return vec15{
		2.0 / 5 * sqrt(numbers::pi) * (l_x + l_y + l_z), // band 0
		0, 0, 4.0 / 7 * sqrt(numbers::pi / 5) * (2 * l_z - (l_x + l_y)), 0, 4.0 / 7 * sqrt(3.0 / 5 * numbers::pi) * (l_x - l_y), // band 2
		0, 0, 0, 0, 2.0 / 35 * sqrt(numbers::pi) * (l_x + l_y + 8.0 / 3 * l_z), 0, -4.0 / 21 * sqrt(numbers::pi / 5) * (l_x - l_y), 0, 2.0 / 3 * sqrt(numbers::pi / 35) * (l_x + l_y) // band 4
	};
}

vec15 odeco_frame(double alpha, double beta, double gamma, double l_x, double l_y, double l_z)
{
	mat15 rotation = rotation_15d(alpha, beta, gamma);
	vec15 ref_frame = reference_frame(l_x, l_y, l_z);
	return rotation * ref_frame;
}

double objective_function(const vec15 y, double alpha, double beta, double gamma, double l_x, double l_y, double l_z)
{
	vec15 f = odeco_frame(alpha, beta, gamma, l_x, l_y, l_z);
	return dot(f - y, f - y);
}

