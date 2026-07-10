
#include "odeco.h"

#include "quaternion.h"

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
	vec15 ref_frame = ref_frame_coords(frame.lambda);
	return rotate(ref_frame, frame.theta);
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
	const double l_x = frame.lambda.x();
	const double l_y = frame.lambda.y();
	const double l_z = frame.lambda.z();
	return vec6{
		1.1102230246251565e-16*l_x*y[7] + l_y*(-3.1381413698186362*y[2] + 1.1298600273745016*y[7] + 1.2811408494623837*y[9]) + l_z*(3.1381413698186362*y[2] + 1.7081877992831784*y[9]),
		l_x*(-1.2811408494623837*y[11] + 1.1298600273745016*y[13] + 3.1381413698186362*y[4]) + 1.1102230246251565e-16*l_y*y[13] - l_z*(1.7081877992831784*y[11] + 3.1381413698186362*y[4]),
		-l_x*(3.1381413698186362*y[1] + 1.5978633742962567*y[6] - 0.60393558820663018*y[8]) - l_y*(-3.1381413698186362*y[1] + 1.5978633742962567*y[6] + 0.60393558820663018*y[8]),
		-C/l_x + 2.7925268031909272*l_x + 0.23935944027350806*l_y + 0.23935944027350806*l_z - 1.4179630807244128*y[0] - 0.20256615438920181*y[10] + 0.30196779410331509*y[12] - 0.39946584357406417*y[14] + 0.90590338230994527*y[3] - 1.5690706849093181*y[5],
		-C/l_y + 0.23935944027350806*l_x + 2.7925268031909272*l_y + 0.23935944027350806*l_z - 1.4179630807244128*y[0] - 0.20256615438920181*y[10] - 0.30196779410331509*y[12] - 0.39946584357406417*y[14] + 0.90590338230994527*y[3] + 1.5690706849093181*y[5],
		-C/l_z + 0.23935944027350806*l_x + 0.23935944027350806*l_y + 2.7925268031909272*l_z - 1.4179630807244128*y[0] - 0.54017641170453823*y[10] - 1.8118067646198905*y[3]
	};
}

mat6 compute_hessian(const vec15& y, const odeco_euler& frame)
{
	const double l_x = frame.lambda.x();
	const double l_y = frame.lambda.y();
	const double l_z = frame.lambda.z();

	return mat6{
	    {2.0*l_x*(2.2204460492503131e-16*y[12] + 5.5511151231257827e-17*y[14]) + 2.0*l_y*(2.0256615438920185*y[10] + 2.4157423528265207*y[12] + 0.79893168714812823*y[14] - 2.7177101469298357*y[3] - 1.5690706849093181*y[5]) + 2.0*l_z*(2.7008820585226916*y[10] + 1.8118067646198912*y[12] + 2.7177101469298357*y[3] + 1.5690706849093181*y[5]),-2.0*l_x*(1.5690706849093181*y[1] + 0.79893168714812823*y[6] - 0.3019677941033152*y[8]) - 1.1102230246251565e-16*l_y*(y[6] - y[8]) + 2.0*l_z*(1.5690706849093181*y[1] + 1.8118067646198905*y[8]),-2.0*l_x*(-0.64057042473119186*y[11] + 0.56493001368725093*y[13] + 1.5690706849093181*y[4]) - 2.0*l_y*(0.64057042473119186*y[11] + 1.6947900410617525*y[13] - 1.5690706849093181*y[4]),0,-3.1381413698186362*y[2] + 1.1298600273745016*y[7] + 1.2811408494623837*y[9],3.1381413698186362*y[2] + 1.7081877992831784*y[9]},
	    {1.1102230246251565e-16*l_x*y[6] + 2.0*l_y*(-1.5690706849093181*y[1] + 0.79893168714812823*y[6] + 0.30196779410331515*y[8]) + 2.0*l_z*(1.5690706849093181*y[1] + 1.8118067646198905*y[8]),2.0*l_x*(2.0256615438920185*y[10] - 2.4157423528265207*y[12] + 0.79893168714812823*y[14] - 2.7177101469298357*y[3] + 1.5690706849093181*y[5]) - 2.0*l_y*(2.2204460492503131e-16*y[12] - 5.5511151231257827e-17*y[14]) + 2.0*l_z*(2.7008820585226916*y[10] - 1.8118067646198912*y[12] + 2.7177101469298357*y[3] - 1.5690706849093181*y[5]),2.0*l_x*(1.5690706849093181*y[2] + 1.6947900410617525*y[7] - 0.64057042473119186*y[9]) + 2.0*l_y*(-1.5690706849093181*y[2] + 0.56493001368725093*y[7] + 0.64057042473119186*y[9]),-1.2811408494623837*y[11] + 1.1298600273745016*y[13] + 3.1381413698186362*y[4],0,-1.7081877992831784*y[11] - 3.1381413698186362*y[4]},
	    {-2.0*l_y*(0.64057042473119186*y[11] + 1.6947900410617522*y[13] - 1.5690706849093181*y[4]) - 2.0*l_z*(0.85409389964158922*y[11] + 1.5690706849093181*y[4]),2.0*l_x*(1.5690706849093181*y[2] + 1.6947900410617522*y[7] - 0.64057042473119186*y[9]) - 2.0*l_z*(1.5690706849093181*y[2] + 0.85409389964158922*y[9]),2.0*l_x*(-0.60393558820663018*y[12] + 3.1957267485925134*y[14] + 3.1381413698186362*y[5]) + 2.0*l_y*(0.60393558820663018*y[12] + 3.1957267485925134*y[14] - 3.1381413698186362*y[5]),-3.1381413698186362*y[1] - 1.5978633742962567*y[6] + 0.60393558820663018*y[8],3.1381413698186362*y[1] - 1.5978633742962567*y[6] - 0.60393558820663018*y[8],0},
	    {0,-1.2811408494623837*y[11] + 1.1298600273745016*y[13] + 3.1381413698186362*y[4],-3.1381413698186362*y[1] - 1.5978633742962567*y[6] + 0.60393558820663018*y[8],C/pow(l_x, 2) + 2.7925268031909272,0.23935944027350806,0.23935944027350806},
	    {-3.1381413698186362*y[2] + 1.1298600273745016*y[7] + 1.2811408494623837*y[9],0,3.1381413698186362*y[1] - 1.5978633742962567*y[6] - 0.60393558820663018*y[8],0.23935944027350806,C/pow(l_y, 2) + 2.7925268031909272,0.23935944027350806},
	    {3.1381413698186362*y[2] + 1.7081877992831784*y[9],-1.7081877992831784*y[11] - 3.1381413698186362*y[4],0,0.23935944027350806,0.23935944027350806,C/pow(l_z, 2) + 2.7925268031909272}
	};

}

odeco_euler closest_seed(const vec15& y, int max_1d_res)
{
	// naive implementation: try out different rotations and pick the best

	odeco_euler best(vec3(0.0, 0.0, 0.0), vec3(1.0, 1.0, 1.0));

	odeco_euler current = best;

	double best_loss = loss(y, best);

	const double step = numbers::pi / (2 * max_1d_res);
	
	for (int i = 0; i < max_1d_res; i++) {
		for (int j = 0; j < max_1d_res; j++) {
			for (int k = 0; k < max_1d_res; k++) {
				current = odeco_euler(step * vec3(i, j, k), current.lambda);
				double current_loss = loss(y, current);
				if (current_loss < best_loss) {
					best_loss = current_loss;
					best = current;
				}
			}
		}
	}

	return best;
}

odeco_mat odeco_frame_project(const vec15& y, int max_1d_res_seed, double tol, int& num_iterations)
{
	odeco_euler current_frame = closest_seed(y, max_1d_res_seed);

	// rotation accumulator
	quaternion q = quaternion::from_euler(current_frame.theta);

	vec15 target = rotate(y, -current_frame.theta);
	current_frame.theta = vec3::Zero();

	int max_iter = 1024;
	num_iterations = max_iter;

	for (int i = 1; i <= max_iter; i++) {
		// compute gradient
		vec6 grad = compute_gradient(target, current_frame);

		// compute hessian
		mat6 hess = compute_hessian(target, current_frame);

		// compute newton step (modify hessian to be pos. def.)
		vec6 newton_step = compute_newton_step(grad, hess);
		
		// compute newton decrement
		double newton_decrement = compute_newton_decrement(grad, newton_step);
		
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
		double t = compute_step_size(target, current_frame, newton_step, -newton_decrement, 1.0, 0.1, 0.2);
		
		// update frame description
		current_frame = update_frame(current_frame, newton_step, t);

		// avoid gimbal locks: absorb rotation into target, keep current_frame.theta = 0
		q = q * quaternion::from_euler(current_frame.theta);
		target = rotate(y, q.conjugate().to_euler());
		current_frame = odeco_euler(current_frame.lambda);
	}
	auto proj = odeco_mat(current_frame.lambda);
	proj.rot = rotation_15d(q.to_euler());
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