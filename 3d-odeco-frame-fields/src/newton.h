#pragma once

#include "linalg.h"

struct odeco_euler;

vec6 compute_newton_step(const vec6 & grad, const mat6 & hess);
double compute_newton_decrement(const vec6& gradient, const vec6& newton_step);
double compute_step_size(const vec15& y, const odeco_euler& frame, const vec6& newton_step, double neg_newton_decr, double t_init, double alpha, double tau);
odeco_euler update_frame(const odeco_euler & frame, const vec6& newton_step, double t);