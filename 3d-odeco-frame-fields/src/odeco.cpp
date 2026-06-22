
#include "odeco.h"

using namespace std;

static const double C = 1e-5;

odeco_euler operator+(const odeco_euler& f0, const odeco_euler& f1)
{
	return odeco_euler(f0.theta + f1.theta, f0.lambda + f1.lambda);
}

odeco_euler operator*(double t, const odeco_euler& f)
{
	return odeco_euler(t * f.theta, t * f.lambda);
}

vec15 ref_frame_coords(const vec3& lambda)
{
	const double l_x = lambda[0], l_y = lambda[1], l_z = lambda[2];

	static const double sqrt_pi = sqrt(numbers::pi);
	static const double sqrt_pi_3 = sqrt(3 * numbers::pi);
	static const double sqrt_5 = sqrt(5);

	return {
		(2 * sqrt_pi / 5) * (l_x + l_y + l_z),
		0,
		0,
		-(4 * sqrt_pi) / (7 * sqrt_5) * (l_x + l_y - 2 * l_z),
		0,
		(4 * sqrt_pi_3) / (7 * sqrt_5) * (l_x - l_y),
		0,
		0,
		0,
		0,
		(2 * sqrt_pi / 105) * (3 * l_x + 3 * l_y + 8 * l_z),
		0,
		(4 * sqrt_pi) / (21 * sqrt_5) * (l_y - l_x),
		0,
		(2 * sqrt_pi) / (3 * sqrt(35)) * (l_x + l_y)
	};
}

vec15 odeco_frame_coords(const odeco_euler& frame)
{
	mat15 rotation = rotation_15d(frame.theta);
	vec15 ref_frame = ref_frame_coords(frame.lambda);
	return block_prod_vec(rotation, ref_frame);
}

vec15 odeco_frame_coords(const odeco_mat& frame)
{
	return block_prod_vec(frame.rot, ref_frame_coords(frame.lambda));
}

double loss(const vec15& y, const odeco_euler& frame)
{
	vec15 f = odeco_frame_coords(frame);
	return (f - y).squaredNorm() - C * (clamped_log(frame.lambda.x()) + clamped_log(frame.lambda.y()) + clamped_log(frame.lambda.z()));
}

double loss(const vec15& y, const odeco_mat& frame)
{
	vec15 f = odeco_frame_coords(frame);
	return (f - y).squaredNorm() - C * (clamped_log(frame.lambda.x()) + clamped_log(frame.lambda.y()) + clamped_log(frame.lambda.z()));
}

vec6 compute_gradient(const vec15& y, const odeco_euler& frame)
{
	vec15 ref_frame = ref_frame_coords(frame.lambda);
	static const mat15 L_x = lie_x_15d();
	static const mat15 L_y = lie_y_15d();
	static const mat15 L_z = lie_z_15d();

	vec6 gradient = vec6::Zero();

	gradient(0) = -2 * y.transpose() * L_x * ref_frame;
	gradient(1) = -2 * y.transpose() * L_y * ref_frame;
	gradient(2) = -2 * y.transpose() * L_z * ref_frame;
	
	gradient.block<3, 1>(3, 0) = 2 * F.transpose() * ref_frame - 2 * F.transpose() * y;

	gradient(3) -= C / frame.lambda.x();
	gradient(4) -= C / frame.lambda.y();
	gradient(5) -= C / frame.lambda.z();

	return gradient;
}

mat6 compute_hessian(const vec15& y, const odeco_euler& frame)
{
	vec15 ref_frame = ref_frame_coords(frame.lambda);
	static const mat15 L_x = lie_x_15d();
	static const mat15 L_y = lie_y_15d();
	static const mat15 L_z = lie_z_15d();
	static const mat15 Lxx = L_x * L_x;
	static const mat15 Lxy = L_x * L_y;
	static const mat15 Lxz = L_x * L_z;
	static const mat15 Lyy = L_y * L_y;
	static const mat15 Lyz = L_y * L_z;
	static const mat15 Lzz = L_z * L_z;

	mat6 H = mat6::Zero();

	const double H_1_12 = -2 * y.transpose() * Lxy * ref_frame;
	const double H_1_13 = -2 * y.transpose() * Lxz * ref_frame;
	const double H_1_23 = -2 * y.transpose() * Lyz * ref_frame;

	mat3 H_1{
		{ -2 * y.transpose() * Lxx * ref_frame, H_1_12,  H_1_13},
		{ H_1_12, -2 * y.transpose() * Lyy * ref_frame,  H_1_23},
		{ H_1_13, H_1_23, -2 * y.transpose() * Lzz * ref_frame },
	};

	mat3 H_2 = mat3::Zero();
	H_2.block<1, 3>(0, 0) = (-2 * F.transpose() * L_x.transpose() * y).transpose();
	H_2.block<1, 3>(1, 0) = (-2 * F.transpose() * L_y.transpose() * y).transpose();
	H_2.block<1, 3>(2, 0) = (-2 * F.transpose() * L_z.transpose() * y).transpose();

	mat3 H_3 = H_2.transpose();

	static const mat3 H_4 = 2 * F.transpose() * F;

	H.block<3, 3>(0, 0) = H_1;
	H.block<3, 3>(0, 3) = H_2;
	H.block<3, 3>(3, 0) = H_3;
	H.block<3, 3>(3, 3) = H_4;

	H(3, 3) += C / (frame.lambda.x() * frame.lambda.x());
	H(4, 4) += C / (frame.lambda.y() * frame.lambda.y());
	H(5, 5) += C / (frame.lambda.z() * frame.lambda.z());

	return H;
}

odeco_euler closest_seed(const vec15& y, int max_1d_res)
{
	// naive implementation: try out different rotations and pick the best

	odeco_euler best(vec3(0.0, 0.0, 0.0), vec3(1.0, 1.0, 1.0));

	odeco_euler current = best;

	double best_loss = loss(y, best);
	double current_loss;

	double step = numbers::pi / (2 * max_1d_res);
	
	for (int i = 0; i < max_1d_res; i++) {
		for (int j = 0; j < max_1d_res; j++) {
			for (int k = 0; k < max_1d_res; k++) {
				current = odeco_euler(step * vec3(i, j, k), current.lambda);
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
odeco_mat odeco_frame_project(const vec15& y, int max_1d_res_seed, double tol, int& num_iterations)
{
	odeco_euler current_frame = closest_seed(y, max_1d_res_seed);

	mat15 rotation = rotation_15d(current_frame.theta);

	vec15 target = rotation.transpose() * y;

	int max_iter = 1024;
	num_iterations = max_iter;

	vec6 grad;
	mat6 hess;
	vec6 newton_step;
	double newton_decrement;
	double t;

	for (int i = 1; i <= max_iter; i++) {
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
			// if the hessian is not pos def then we are in a saddle point

			Eigen::LLT<mat6> decomp = llt(hess);

			if (decomp.info() == Eigen::Success) {
				//std::cout << "break standard" << std::endl;
				num_iterations = i;
				break;
			}
			//std::cout << "saddle point detected" << std::endl;
			// hess is not pos def.
			// find smallest eigenvalue and corresponding eigenvector
			// the newton step will be 0.5 * this eigenvector
			vec6 eivec = eivector_smallest_ev(hess);
			newton_step = 0.5 * eivec;

			newton_decrement = compute_newton_decrement(grad, newton_step);
			if (newton_decrement / 2 <= tol) {
				//std::cout << "break in saddle point" << std::endl;
				num_iterations = i;
				break;
			}
		}
		
		// line search step size
		t = compute_step_size(target, current_frame, newton_step, -newton_decrement, 1.0, 0.1, 0.5);
		
		// update frame description
		current_frame = update_frame(current_frame, newton_step, t);

		// avoid gimbal locks: absorb rotation into target, keep current_frame.theta = 0
		rotation = block_prod_mat(rotation, rotation_15d(current_frame.theta));
		target = block_prod_vec(rotation.transpose(), y);
		current_frame = odeco_euler(current_frame.lambda);
	}
	odeco_mat proj = odeco_mat(current_frame.lambda);
	proj.rot = rotation;
	return proj;
}

odeco_euler::odeco_euler(const vec3& theta, const vec3& lambda)
{
	this->theta = theta;
	this->lambda = lambda;
}

odeco_euler::odeco_euler(const vec3& lambda) : odeco_euler(vec3(0, 0, 0), lambda)
{
}

odeco_mat::odeco_mat(const vec3& theta, const vec3& lambda)
{
	this->rot = rotation_15d(theta);
	this->lambda = lambda;
}

odeco_mat::odeco_mat(const vec3& lambda) : odeco_mat(vec3::Zero(), lambda)
{
}
