#include "odeco_aligned.h"

using namespace std;

static const double C = 1e-5;

odeco_euler_aligned::odeco_euler_aligned(double theta, const vec2& lambda)
{
	this->theta = theta;
	this->lambda = lambda;
}

odeco_euler_aligned::odeco_euler_aligned(const vec2& lambda) : odeco_euler_aligned(0, lambda)
{
}

odeco_quat_aligned::odeco_quat_aligned(const quaternion &rot, const vec2 &lambda) : rot(rot), lambda(lambda)
{
}

odeco_quat_aligned::odeco_quat_aligned(double theta, const vec2& lambda) : rot(quaternion::rotate_z(theta)), lambda(lambda)
{
}

odeco_quat_aligned::odeco_quat_aligned(const vec2& lambda) : odeco_quat_aligned(0, lambda)
{
}

odeco_euler_aligned operator+(const odeco_euler_aligned& f0, const odeco_euler_aligned& f1)
{
	return odeco_euler_aligned(f0.theta + f1.theta, f0.lambda + f1.lambda);
}

odeco_euler_aligned operator*(double t, const odeco_euler_aligned& f)
{
	return odeco_euler_aligned(t * f.theta, t * f.lambda);
}

vec15 ref_frame_coords_aligned(const vec2& lambda)
{
    vec3 l = vec3(lambda.x(), lambda.y(), 1.0);
    return ref_frame_coords(l);
}

vec15 odeco_frame_coords(const odeco_euler_aligned& frame)
{
    mat15 rotation = rotation_z_15d(frame.theta);
    vec15 ref_frame = ref_frame_coords_aligned(frame.lambda);
    return block_prod_vec(rotation, ref_frame);
}

vec15 odeco_frame_coords(const odeco_quat_aligned& frame)
{
	return rotate_z(ref_frame_coords_aligned(frame.lambda), frame.rot.to_euler().z());
}

double loss(const vec15& y, const odeco_euler_aligned& frame)
{
    vec15 f = odeco_frame_coords(frame);
    return (f - y).squaredNorm() - C * ( clamped_log(frame.lambda.x()) + clamped_log(frame.lambda.y()));
}

double loss(const vec15& y, const odeco_quat_aligned& frame)
{
    vec15 f = odeco_frame_coords(frame);
    return (f - y).squaredNorm() - C * ( clamped_log(frame.lambda.x()) + clamped_log(frame.lambda.y()));
}

vec3 compute_gradient(const vec15& y, const odeco_euler_aligned& frame)
{
	const auto l_x = frame.lambda.x();
	const auto l_y = frame.lambda.y();

	return vec3{
		-l_x*(3.1381413698186362*y[1] + 1.5978633742962567*y[6] - 0.60393558820663018*y[8]) - l_y*(-3.1381413698186362*y[1] + 1.5978633742962567*y[6] + 0.60393558820663018*y[8]),
		-C*l_x + 2.7925268031909272*l_x + 0.23935944027350806*l_y - 1.4179630807244128*y[0] - 0.20256615438920181*y[10] + 0.30196779410331509*y[12] - 0.39946584357406417*y[14] + 0.90590338230994527*y[3] - 1.5690706849093181*y[5] + 0.23935944027350806,
		-C*l_y + 0.23935944027350806*l_x + 2.7925268031909272*l_y - 1.4179630807244128*y[0] - 0.20256615438920181*y[10] - 0.30196779410331509*y[12] - 0.39946584357406417*y[14] + 0.90590338230994527*y[3] + 1.5690706849093181*y[5] + 0.23935944027350806
	};
}

mat3 compute_hessian(const vec15& y, const odeco_euler_aligned& frame)
{
	const auto l_x = frame.lambda.x();
	const auto l_y = frame.lambda.y();

	return mat3{
	    {l_x*(-1.2078711764132604*y[12] + 6.3914534971850268*y[14] + 6.2762827396372725*y[5]) + l_y*(1.2078711764132604*y[12] + 6.3914534971850268*y[14] - 6.2762827396372725*y[5]),-3.1381413698186362*y[1] - 1.5978633742962567*y[6] + 0.60393558820663018*y[8],3.1381413698186362*y[1] - 1.5978633742962567*y[6] - 0.60393558820663018*y[8]},
		{-3.1381413698186362*y[1] - 1.5978633742962567*y[6] + 0.60393558820663018*y[8],C/pow(l_x, 2) + 2.7925268031909272,0.23935944027350806},
		{3.1381413698186362*y[1] - 1.5978633742962567*y[6] - 0.60393558820663018*y[8],0.23935944027350806,C/pow(l_y, 2) + 2.7925268031909272}
	};
}

odeco_quat_aligned odeco_frame_project_aligned(const vec15& y, const vec3& d, double tol, int& num_iterations)
{
	if (d.cwiseAbs().maxCoeff() < 1e-2) {
		cerr << "Alignement axis is the zero vector - please provide a different alignement vector" << endl;
		std::exit(1);
	}

	const auto z_to_d = quaternion::z_to_dir(d);
	const auto target = rotate(y, z_to_d.conjugate().to_euler());

	const auto z_aligned_target = vec15(
		0.69999999999999996 * target(0) + 0.10000000000000001 * target(10) - 0.44721359549995798 * target(3) + 0.59081795030183859,
		0.9642857142857143 * target(1) - 0.18557687223952254 * target(8),
		0,
		-0.44721359549995798 * target(0) - 0.063887656499993992 * target(10) + 0.28571428571428575 * target(3) + 0.98139533083577413,
		0,
		-0.18557687223952254 * target(12) + 0.9642857142857143 * target(5),
		target(6),
		0,
		-0.18557687223952254 * target(1) + 0.035714285714285712 * target(8),
		0,
		0.10000000000000001 * target(0) + 0.014285714285714285 * target(10) - 0.063887656499993992 * target(3) + 0.25320769298650225,
		0,
		0.035714285714285712 * target(12) - 0.18557687223952254 * target(5),
		0,
		target(14)
	);

	odeco_quat_aligned z_aligned_proj = odeco_frame_project_aligned_newton(z_aligned_target, tol, num_iterations);
	z_aligned_proj.rot = z_to_d * z_aligned_proj.rot;
	return z_aligned_proj;

}

odeco_quat_aligned odeco_frame_project_aligned_newton(const vec15& y, double tol, int& num_iterations)
{
	odeco_euler_aligned current_frame = odeco_euler_aligned(0, vec2(1, 1));

	mat15 rot = rotation_z_15d(current_frame.theta);
	auto q = quaternion::rotate_z(current_frame.theta);

	auto target = rotate(y, q.conjugate().to_euler());

	int max_iter = 1024;
	num_iterations = max_iter;

	for (int i = 1; i <= max_iter; i++) {
		// compute gradient
		vec3 grad = compute_gradient(target, current_frame);

		// compute hessian
		mat3 hess = compute_hessian(target, current_frame);

		// compute newton step
		vec3 newton_step = compute_newton_step(grad, hess);

		// compute newton decrement
		double newton_decrement = compute_newton_decrement(grad, newton_step);

		// stop if dec < tol
		if (newton_decrement / 2 <= tol) {
			Eigen::LLT<mat3> decomp = llt(hess);

			if (decomp.info() == Eigen::Success) {
				num_iterations = i;
				break;
			}

			vec3 eivec = eivector_smallest_ev(hess);
			newton_step = 0.5 * eivec;

			newton_decrement = compute_newton_decrement(grad, newton_step);
			if (newton_decrement / 2 <= tol) {
				num_iterations = i;
				break;
			}
		}

		// line search step size
		double t = compute_step_size(target, current_frame, newton_step, -newton_decrement, 1.0, 0.1, 0.2);

		// update frame description
		current_frame = update_frame(current_frame, newton_step, t);

		// avoid gimbal lock
		q = q * quaternion::rotate_z(current_frame.theta);
		target = rotate_z(y, q.conjugate().to_euler().z());
		current_frame = odeco_euler_aligned(current_frame.lambda);

	}
	return odeco_quat_aligned(q, current_frame.lambda);
}