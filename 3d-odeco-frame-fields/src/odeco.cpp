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

odeco_frame_description operator+(const odeco_frame_description& f0, const odeco_frame_description& f1)
{
	return odeco_frame_description{
		.theta = f0.theta + f1.theta,
		.lambda = f0.lambda + f1.lambda
	};
}

odeco_frame_description operator*(double t, const odeco_frame_description& f)
{
	return odeco_frame_description{
		.theta = t * f.theta,
		.lambda = t * f.lambda
	};
}

vec15 reference_frame(const vec3& lambda)
{
	return F * lambda;
}

vec15 odeco_frame(const odeco_frame_description& frame)
{
	mat15 rotation = rotation_15d(frame.theta);
	vec15 ref_frame = reference_frame(frame.lambda);
	return rotation * ref_frame;
}

double loss(const vec15& y, const odeco_frame_description& frame)
{
	vec15 f = odeco_frame(frame);
	return dot(f - y, f - y);
}

vec6 compute_gradient(const vec15& y, const odeco_frame_description& frame)
{
	vec15 ref_frame = reference_frame(frame.lambda);
	mat15 L_x = lie_x_15d();
	mat15 L_y = lie_y_15d();
	mat15 L_z = lie_z_15d();
	mat15 R_x = rotation_x_15d(frame.theta.x());
	mat15 R_y = rotation_y_15d(frame.theta.y());
	mat15 R_z = rotation_z_15d(frame.theta.z());

	vec6 gradient = vec6::Zero();

	gradient(0) = -2 * y.transpose() * L_x * R_x * R_y * R_z * ref_frame;
	gradient(1) = -2 * y.transpose() * R_x * L_y * R_y * R_z * ref_frame;
	gradient(2) = -2 * y.transpose() * R_x * R_y * L_z * R_z * ref_frame;
	
	gradient.block<3, 1>(3, 0) = 2 * F.transpose() * ref_frame - 2 * F.transpose() * (R_x * R_y * R_z).transpose() * y;

	return gradient;
}

mat6 compute_hessian(const vec15& y, const odeco_frame_description& frame)
{
	vec15 ref_frame = reference_frame(frame.lambda);
	mat15 L_x = lie_x_15d();
	mat15 L_y = lie_y_15d();
	mat15 L_z = lie_z_15d();
	mat15 R_x = rotation_x_15d(frame.theta.x());
	mat15 R_y = rotation_y_15d(frame.theta.y());
	mat15 R_z = rotation_z_15d(frame.theta.z());

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

odeco_frame_description closest_seed(const vec15& y, int max_1d_res)
{
	// naive implementation: try out different rotations and pick the best one

	odeco_frame_description best = {
		.theta = vec3(0.0, 0.0, 0.0),
		.lambda = vec3(1.0, 1.0, 1.0)
	};

	odeco_frame_description current = best;

	double best_loss = loss(y, best);
	double current_loss;

	double step = numbers::pi / max_1d_res;
	
	for (int i = 0; i < max_1d_res; i++) {
		for (int j = 0; j < max_1d_res; j++) {
			for (int k = 0; k < max_1d_res; k++) {
				current.theta = (step / numbers::pi) * vec3(i, j, k);
				current_loss = loss(y, current);
				if (current_loss < best_loss) {
					best_loss = current_loss;
					best = current;
				}
			}
		}
	}

	return best;
}

odeco_frame_description odeco_frame_project(const vec15& y, int max_1d_res_seed, double tol)
{
	odeco_frame_description current_frame_description = closest_seed(y, max_1d_res_seed);

	vec6 current_x = vec6{
		current_frame_description.theta(0),
		current_frame_description.theta(1),
		current_frame_description.theta(2),
		current_frame_description.lambda(0),
		current_frame_description.lambda(1),
		current_frame_description.lambda(2),
	};

	int max_iter = 1024;

	vec15 current_frame;
	vec6 grad;
	mat6 hess;
	vec6 newton_step;
	double newton_decrement;
	double t;

	for (int i = 0; i < max_iter; i++) {

		// compute current frame
		current_frame = odeco_frame(current_frame_description);

		// compute gradient
		grad = compute_gradient(y, current_frame_description);

		// compute hessian
		hess = compute_hessian(y, current_frame_description);

		// compute newton step (modify hessian to be pos. def.)
		newton_step = compute_newton_step(grad, hess);
		
		// compute newton decrement
		newton_decrement = compute_newton_decrement(grad, newton_step);
		
		// stop if dec < tol
		if (newton_decrement / 2 <= tol) break;
		
		// line search step size
		t = compute_step_size(y, current_frame_description, newton_step, grad, 1.0, 0.1, 0.9);
		
		// update frame description
		current_frame_description = update_frame_description(current_frame_description, newton_step, t);
	}

	return current_frame_description;
}
