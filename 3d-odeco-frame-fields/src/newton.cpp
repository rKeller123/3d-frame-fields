#include "newton.h"
#include "odeco.h"

vec6 compute_newton_step(const vec6& grad, const mat6& hess)
{
	double epsilon = 0.1;

	vec6 eigenvalues = eivals(hess);
	double min_ev = eigenvalues.minCoeff();

	min_ev = std::min(min_ev - epsilon, 0.0);

	mat6 H_mod = hess - min_ev * mat6::Identity();

	// TODO: use cholesky decomposition and solve for the inverse
	return -H_mod.inverse() * grad;
}

double compute_newton_decrement(const vec6& gradient, const vec6& newton_step)
{
	return -dot(gradient, newton_step);
}

double compute_step_size(const vec15& y, const odeco_frame& frame, const vec6& dx, const vec6& grad, double t_init, double alpha, double tau)
{
	odeco_frame dx_frame{
		.theta = dx.block<3, 1>(0, 0),
		.lambda = dx.block<3, 1>(3, 0)
	};

	double t = t_init;

	while (loss(y, frame + t * dx_frame) > loss(y, frame) + t * alpha * dot(grad, dx)) {
		t *= tau;
	}

	return t;
}

odeco_frame update_frame_description(const odeco_frame& frame, const vec6& newton_step, double t)
{
	odeco_frame newton_step_frame{
		.theta = newton_step.block<3, 1>(0, 0),
		.lambda = newton_step.block<3, 1>(3, 0)
	};
	return frame + t * newton_step_frame;
}
