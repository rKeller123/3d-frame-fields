#include "odeco.h"


using namespace std;

mat15x3 F{
		{sqrt(numbers::pi) * 2.0 / 5, sqrt(numbers::pi) * 2.0 / 5, sqrt(numbers::pi) * 2.0 / 5},

		{0, 0, 0},
		{0, 0, 0},
		{-(4.0 * sqrt(numbers::pi)) / (7 * sqrt(5)), -(4.0 * sqrt(numbers::pi)) / (7 * sqrt(5)), 2 * (4.0 * sqrt(numbers::pi)) / (7 * sqrt(5))},
		{0, 0, 0},
		{(4 * sqrt(3 * numbers::pi)) / (7 * sqrt(5)), -(4 * sqrt(3 * numbers::pi)) / (7 * sqrt(5)), 0},

		{0, 0, 0},
		{0, 0, 0},
		{0, 0, 0},
		{0, 0, 0},
		{2 * sqrt(numbers::pi) / 35, 2 * sqrt(numbers::pi) / 35, 16 * sqrt(numbers::pi) / 105},
		{0, 0, 0},
		{-(4 * sqrt(numbers::pi)) / (21 * sqrt(5)), (4 * sqrt(numbers::pi)) / (21 * sqrt(5)), 0},
		{0, 0, 0},
		{(2 * sqrt(numbers::pi)) / (3 * sqrt(35)), (2 * sqrt(numbers::pi)) / (3 * sqrt(35)), 0}
};

vec15 reference_frame(double l_x, double l_y, double l_z)
{

	vec3 lambda{l_x, l_y, l_z};
	return F * lambda;
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

vec6 gradient(const vec15 y, double alpha, double beta, double gamma, double l_x, double l_y, double l_z)
{
	vec15 ref_frame = reference_frame(l_x, l_y, l_z);
	mat15 L_x = lie_x_15d();
	mat15 L_y = lie_y_15d();
	mat15 L_z = lie_z_15d();
	mat15 R_x = rotation_x_15d(alpha);
	mat15 R_y = rotation_y_15d(beta);
	mat15 R_z = rotation_z_15d(gamma);

	vec6 gradient = vec6::Zero();

	gradient(0) = -2 * y.transpose() * L_x * R_x * R_y * R_z * ref_frame;
	gradient(1) = -2 * y.transpose() * R_x * L_y * R_y * R_z * ref_frame;
	gradient(2) = -2 * y.transpose() * R_x * R_y * L_z * R_z * ref_frame;
	
	gradient.block<3, 1>(3, 0) = 2 * F.transpose() * ref_frame - 2 * F.transpose() * (R_x * R_y * R_z).transpose() * y;

	return gradient;
}

mat6 hessian(const vec15 y, double alpha, double beta, double gamma, double l_x, double l_y, double l_z)
{
	vec15 ref_frame = reference_frame(l_x, l_y, l_z);
	mat15 L_x = lie_x_15d();
	mat15 L_y = lie_y_15d();
	mat15 L_z = lie_z_15d();
	mat15 R_x = rotation_x_15d(alpha);
	mat15 R_y = rotation_y_15d(beta);
	mat15 R_z = rotation_z_15d(gamma);

	mat6 H = mat6::Zero();

	mat3 H_1{
		{-2 * y.transpose() * L_x * L_x * R_x * R_y * R_z * ref_frame, -2 * y.transpose() * L_x * R_x * L_y * R_y * R_z * ref_frame, -2 * y.transpose() * L_x * R_x * R_y * L_z * R_z * ref_frame},
		{-2 * y.transpose() * L_x * R_x * L_y * R_y * R_z * ref_frame, -2 * y.transpose() * R_x * L_y * L_y * R_y * R_z * ref_frame, -2 * y.transpose() * R_x * L_y * R_y * L_z * R_z * ref_frame},
		{-2 * y.transpose() * L_x * R_x * R_y * L_z * R_z * ref_frame, -2 * y.transpose() * R_x * L_y * R_y * L_z * R_z * ref_frame, -2 * y.transpose() * R_x * R_y * L_z * L_z * R_z * ref_frame},
	};

	mat3 H_2 = mat3::Zero();
	H_2.block<1, 3>(0, 0) = (-2 * F.transpose() * (L_x * R_x * R_y * R_z).transpose() * y).transpose();
	H_2.block<1, 3>(1, 0) = (-2 * F.transpose() * (R_x * L_y * R_y * R_z).transpose() * y).transpose();
	H_2.block<1, 3>(2, 0) = (-2 * F.transpose() * (R_x * R_y * L_z * R_z).transpose() * y).transpose();

	mat3 H_3 = H_2.transpose();

	mat3 H_4 = 2 * F.transpose() * F;

	H.block<3, 3>(0, 0) = H_1;
	H.block<3, 3>(0, 3) = H_2;
	H.block<3, 3>(3, 0) = H_3;
	H.block<3, 3>(3, 3) = H_4;

	return H;
}

