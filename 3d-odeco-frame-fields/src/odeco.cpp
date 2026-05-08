
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

odeco_frame operator+(const odeco_frame& f0, const odeco_frame& f1)
{
	return odeco_frame{
		.theta = f0.theta + f1.theta,
		.lambda = f0.lambda + f1.lambda
	};
}

odeco_frame operator*(double t, const odeco_frame& f)
{
	return odeco_frame{
		.theta = t * f.theta,
		.lambda = t * f.lambda
	};
}

vec15 ref_frame_coords(const vec3& lambda)
{
	return F * lambda;
}

vec15 odeco_frame_coords(const odeco_frame& frame)
{
	mat15 rotation = rotation_15d(frame.theta);
	vec15 ref_frame = ref_frame_coords(frame.lambda);
	return rotation * ref_frame;
}

double loss(const vec15& y, const odeco_frame& frame)
{
	vec15 f = odeco_frame_coords(frame);
	return loss(y, f);
}

double loss(const vec15& y, const vec15& frame_coords)
{
	return dot(frame_coords - y, frame_coords - y);
}

vec6 compute_gradient(const vec15& y, const odeco_frame& frame)
{
	vec15 ref_frame = ref_frame_coords(frame.lambda);
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

mat6 compute_hessian(const vec15& y, const odeco_frame& frame)
{
	vec15 ref_frame = ref_frame_coords(frame.lambda);
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

odeco_frame closest_seed(const vec15& y, int max_1d_res)
{
	// naive implementation: try out different rotations and pick the best

	odeco_frame best = {
		.theta = vec3(0.0, 0.0, 0.0),
		.lambda = vec3(1.0, 1.0, 1.0)
	};

	odeco_frame current = best;

	double best_loss = loss(y, best);
	double current_loss;

	double step = numbers::pi / (2 * max_1d_res);
	
	for (int i = 0; i < max_1d_res; i++) {
		for (int j = 0; j < max_1d_res; j++) {
			for (int k = 0; k < max_1d_res; k++) {
				current.theta = step * vec3(i, j, k);
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

// TODO: consider non-negativity constraints for lambda
vec15 odeco_frame_project(const vec15& y, int max_1d_res_seed, double tol, int& num_iterations)
{
	odeco_frame seed = closest_seed(y, max_1d_res_seed);

	mat15 rotation = rotation_15d(seed.theta);

	vec15 target = rotation.transpose() * y;
	odeco_frame current_frame{
		.theta = vec3(0, 0, 0),
		.lambda = vec3(1, 1, 1)
	};

	int max_iter = 1024;
	num_iterations = max_iter;

	vec15 current_frame_coords;
	vec6 grad;
	mat6 hess;
	vec6 newton_step;
	mat15 rotation_step;
	double newton_decrement;
	double t;

	for (int i = 0; i < max_iter; i++) {
		// compute gradient
		grad = compute_gradient(target, current_frame);

		// compute hessian
		hess = compute_hessian(target, current_frame);

		// compute newton step (modify hessian to be pos. def.)
		newton_step = compute_newton_step(grad, hess);
		
		// compute newton decrement
		newton_decrement = compute_newton_decrement(grad, newton_step);
		
		// stop if dec < tol
		if (newton_decrement / 2 <= tol) {
			num_iterations = i;
			break;
		}
		
		// line search step size
		t = compute_step_size(target, current_frame, newton_step, grad, 1.0, 0.1, 0.9);
		
		// update frame description
		current_frame = update_frame(current_frame, newton_step, t);

		// avoid gimbal locks: absorb rotation into target, keep current_frame.theta = 0
		rotation_step = rotation_15d(current_frame.theta);
		target = rotation_step.transpose() * target;
		rotation = rotation * rotation_step;
		current_frame.theta = vec3(0, 0, 0);
	}

	return rotation * odeco_frame_coords(current_frame);
}
