#pragma once

#include "linalg.h"
#include "newton.h"
#include <numbers>
#include <iostream>

struct odeco_euler {
	vec3 theta;
	vec3 lambda;
	odeco_euler(const vec3& theta, const vec3& lambda);
	odeco_euler(const vec3& lambda);
};

struct odeco_mat {
	mat15 rot;
	vec3 lambda;

	odeco_mat(const vec3& theta, const vec3& lambda);
	odeco_mat(const vec3& lambda);
};

odeco_euler operator+(const odeco_euler& f0, const odeco_euler& f1);
odeco_euler operator*(double t, const odeco_euler& f);

vec15 ref_frame_coords(const vec3& lambda);
vec15 odeco_frame_coords(const odeco_euler& frame);
vec15 odeco_frame_coords(const odeco_mat& frame);

double loss(const vec15& y, const odeco_euler& frame);
double loss(const vec15& y, const odeco_mat& frame);

vec6 compute_gradient(const vec15& y, const odeco_euler& frame);

mat6 compute_hessian(const vec15& y, const odeco_euler& frame);

odeco_euler closest_seed(const vec15& y, int max_1d_res);

odeco_mat odeco_frame_project(const vec15& y, int max_1d_res_seed, double tol, int& num_iterations);
odeco_mat odeco_frame_project_aligned(const vec15& y, const vec3& d, int max_1d_res_seed, double tol, int& num_iterations);