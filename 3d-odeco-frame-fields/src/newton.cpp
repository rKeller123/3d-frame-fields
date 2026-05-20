#include "newton.h"
#include "odeco.h"

vec6 compute_newton_step(const vec6& grad, const mat6& hess)
{
	double epsilon = 0.1;
	double gamma = 10.0;
	double delta = std::max(abs(0.000015 * hess.trace()), epsilon);
	double m = delta;

	mat6 H_mod = hess;

	Eigen::LLT<mat6> decomp = llt(H_mod);

	while (decomp.info() != Eigen::Success) {
		H_mod = H_mod + m * mat6::Identity();
		m = gamma * m;
		decomp = llt(H_mod);
	}

	mat6 L = decomp.matrixL();
	vec6 y = L.triangularView<Eigen::Lower>().solve(-grad);
	return L.transpose().triangularView<Eigen::Upper>().solve(y);
}

double compute_newton_decrement(const vec6& gradient, const vec6& newton_step)
{
	return -dot(gradient, newton_step);
}

double compute_step_size(const vec15& y, const odeco_euler& frame, const vec6& newton_step, double neg_newton_decr, double t_init, double alpha, double tau)
{
	odeco_euler dx_frame(newton_step.block<3, 1>(0, 0), newton_step.block<3, 1>(3, 0));

	double loss_y_frame = loss(y, frame);
	double alpha_dot = alpha * neg_newton_decr;

	double t = t_init;

	while (loss(y, frame + t * dx_frame) > loss_y_frame + t * alpha_dot) {
		t *= tau;
	}

	return t;
}

odeco_euler update_frame(const odeco_euler& frame, const vec6& newton_step, double t)
{
	odeco_euler newton_step_frame(newton_step.block<3, 1>(0, 0), newton_step.block<3, 1>(3, 0));
	odeco_euler update = frame + t * newton_step_frame;
	return update;
}
