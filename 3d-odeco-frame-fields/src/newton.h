#pragma once

#include "linalg.h"

struct odeco_euler;
struct odeco_euler_aligned;

template<typename vec, typename mat>
vec compute_newton_step(const vec & grad, const mat & hess)
{
    double epsilon = 0.1;
    double gamma = 10.0;
    double delta = std::max(abs(0.000015 * hess.trace()), epsilon);
    double m = delta;

    mat H_mod = hess;

    Eigen::LLT<mat> decomp = llt(H_mod);

    int n_fails = 0;

    while (decomp.info() != Eigen::Success) {
        H_mod = H_mod + m * mat::Identity();
        m = gamma * m;
        decomp = llt(H_mod);
        n_fails++;
    }

    //std::cout << "N failed cholesky: " << n_fails << std::endl;

    return decomp.solve(-grad);
}

template<typename vec>
double compute_newton_decrement(const vec& gradient, const vec& newton_step)
{
    return -dot(gradient, newton_step);
}

double compute_step_size(const vec15& y, const odeco_euler& frame, const vec6& newton_step, double neg_newton_decr, double t_init, double alpha, double tau);
double compute_step_size(const vec15& y, const odeco_euler_aligned& frame, const vec3& newton_step, double neg_newton_decr, double t_init, double alpha, double tau);

odeco_euler update_frame(const odeco_euler& frame, const vec6& newton_step, double t);
odeco_euler_aligned update_frame(const odeco_euler_aligned& frame, const vec3& newton_step, double t);