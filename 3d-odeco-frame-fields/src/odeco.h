#pragma once

#include "linalg.h"
#include "newton.h"
#include <numbers>

struct odeco_frame_description {
	vec3 theta;
	vec3 lambda;
};

odeco_frame_description operator+(const odeco_frame_description& f0, const odeco_frame_description& f1);
odeco_frame_description operator*(double t, const odeco_frame_description& f);

vec15 reference_frame(const vec3& lambda);
vec15 odeco_frame(const odeco_frame_description& frame);

double loss(const vec15& y, const odeco_frame_description& frame);

vec6 compute_gradient(const vec15& y, const odeco_frame_description& frame);

mat6 compute_hessian(const vec15& y, const odeco_frame_description& frame);

odeco_frame_description closest_seed(const vec15& y, int max_1d_res);

odeco_frame_description odeco_frame_project(const vec15& y, int max_1d_res_seed, double tol);