#pragma once

#include "linalg.h"

struct odeco_frame;

vec6 compute_newton_step(const vec6 & grad, const mat6 & hess);
double compute_newton_decrement(const vec6& gradient, const vec6& newton_step);
double compute_step_size(const vec15& y, const odeco_frame& frame, const vec6& dx, const vec6& grad, double t_init, double alpha, double tau);
odeco_frame update_frame_description(const odeco_frame& frame, const vec6& newton_step, double t);