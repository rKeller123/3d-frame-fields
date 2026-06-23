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

odeco_mat_aligned::odeco_mat_aligned(double theta, const vec2& lambda)
{
	this->rot = rotation_z_15d(theta);
	this->lambda = lambda;
}

odeco_mat_aligned::odeco_mat_aligned(const vec2& lambda) : odeco_mat_aligned(0, lambda)
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

vec15 odeco_frame_coords(const odeco_mat_aligned& frame)
{
    return block_prod_vec(frame.rot, ref_frame_coords_aligned(frame.lambda));
}

double loss(const vec15& y, const odeco_euler_aligned& frame)
{
    vec15 f = odeco_frame_coords(frame);
    return (f - y).squaredNorm() - C * ( clamped_log(frame.lambda.x()) + clamped_log(frame.lambda.y()));
}

double loss(const vec15& y, const odeco_mat_aligned& frame)
{
    vec15 f = odeco_frame_coords(frame);
    return (f - y).squaredNorm() - C * ( clamped_log(frame.lambda.x()) + clamped_log(frame.lambda.y()));
}

vec3 compute_gradient(const vec15& y, const odeco_euler_aligned& frame)
{
    vec15 ref_frame = ref_frame_coords_aligned(frame.lambda);
    static const mat15 L_z = lie_z_15d();

    vec3 gradient = vec3::Zero();

    gradient(0) = -2 * y.transpose() * L_z * ref_frame;
    
    vec3 grad_lambda_full = 2 * F.transpose() * ref_frame - 2 * F.transpose() * y;

    gradient(1) = grad_lambda_full.x() - C / frame.lambda.x();
    gradient(2) = grad_lambda_full.y() - C / frame.lambda.y();

    return gradient;
}

mat3 compute_hessian(const vec15& y, const odeco_euler_aligned& frame)
{
    vec15 ref_frame = ref_frame_coords_aligned(frame.lambda);
	static const mat15 L_z = lie_z_15d();
	static const mat15 Lzz = L_z * L_z;

    mat3 H = mat3::Zero();

    H(0, 0) = -2 * y.transpose() * Lzz * ref_frame;

    vec3 hess_theta_full = -2 * F.transpose() * L_z.transpose() * y;

    H(0, 1) = hess_theta_full.x();
    H(0, 2) = hess_theta_full.y();
    H(1, 0) = hess_theta_full.x();
    H(2, 0) = hess_theta_full.y();

    mat3 hess_lambda_full = 2 * F.transpose() * F;

    H(1, 1) = hess_lambda_full(0, 0) + C / (frame.lambda.x() * frame.lambda.x());
    H(1, 2) = hess_lambda_full(0, 1);
    H(2, 1) = hess_lambda_full(1, 0);
    H(2, 2) = hess_lambda_full(1, 1) + C / (frame.lambda.y() * frame.lambda.y());

    return H;
}

odeco_mat_aligned odeco_frame_project_aligned(const vec15& y, const vec3& d, double tol, int& num_iterations)
{
	if (d.cwiseAbs().maxCoeff() < 1e-2) {
		cerr << "Alignement axis is the zero vector - please provide a different alignement vector" << endl;
		std::exit(1);
	}

	vec3 d_normalized = d.normalized();

	vec3 z = vec3(0, 0, 1);

	vec15 target = y;
	mat15 R_z_to_d;
	bool isZeqD = (z - d_normalized).cwiseAbs().maxCoeff() < 1e-2;
	if (!isZeqD) {
		double angle = acos(dot(z, d_normalized));
		vec3 axis = cross(z, d_normalized);
		R_z_to_d = rotation_from_axis_angle(axis, angle);
		target = R_z_to_d.transpose() * y;

	}

	vec15 z_aligned_target = vec15(
		0.69999999999999996 * target(0) + 0.10000000000000001 * target(10) - 0.44721359549995798 * target(3) + 0.59081795030183859,
		0.9642857142857143 * target(1) - 0.18557687223952254 * y(8),
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

	odeco_mat_aligned z_aligned_proj = odeco_frame_project_aligned_newton(z_aligned_target, tol, num_iterations);
	
	if (!isZeqD) {
		z_aligned_proj.rot = block_prod_mat(R_z_to_d, z_aligned_proj.rot);
	}
	return z_aligned_proj;

}

odeco_mat_aligned odeco_frame_project_aligned_newton(const vec15& y, double tol, int& num_iterations)
{

    odeco_euler_aligned current_frame = odeco_euler_aligned(0, vec2(1, 1));

    mat15 rotation = rotation_z_15d(current_frame.theta);

    vec15 target = rotation.transpose() * y;

	int max_iter = 1024;
    num_iterations = max_iter;

    vec3 grad;
    mat3 hess;
    vec3 newton_step;
    double newton_decrement;
    double t;

	for (int i = 1; i <= max_iter; i++) {
		// compute gradient
        grad = compute_gradient(target, current_frame);

		// compute hessian
		hess = compute_hessian(target, current_frame);
        
        // compute newton step
        newton_step = compute_newton_step(grad, hess);

		// compute newton decrement
        newton_decrement = compute_newton_decrement(grad, newton_step);

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
        t = compute_step_size(target, current_frame, newton_step, -newton_decrement, 1.0, 0.1, 0.5);
        
        // update frame description
        current_frame = update_frame(current_frame, newton_step, t);

        // avoid gimbal lock
        rotation = block_prod_mat(rotation, rotation_z_15d(current_frame.theta));
        target = block_prod_vec(rotation.transpose(), y);
        current_frame = odeco_euler_aligned(current_frame.lambda);

	}

    odeco_mat_aligned proj = odeco_mat_aligned(current_frame.lambda);
    proj.rot = rotation;

    return proj;
}