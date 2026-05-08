#pragma once

#include "linalg.h"
#include "newton.h"
#include <numbers>
#include <iostream>

struct odeco_frame {
	vec3 theta;
	vec3 lambda;
};

odeco_frame operator+(const odeco_frame& f0, const odeco_frame& f1);
odeco_frame operator*(double t, const odeco_frame& f);

vec15 ref_frame_coords(const vec3& lambda);
vec15 odeco_frame_coords(const odeco_frame& frame);

double loss(const vec15& y, const odeco_frame& frame);
double loss(const vec15& y, const vec15& frame_coords);

vec6 compute_gradient(const vec15& y, const odeco_frame& frame);

mat6 compute_hessian(const vec15& y, const odeco_frame& frame);

odeco_frame closest_seed(const vec15& y, int max_1d_res);

vec15 odeco_frame_project(const vec15& y, int max_1d_res_seed, double tol, int& num_iterations);