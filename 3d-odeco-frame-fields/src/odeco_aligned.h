#pragma once

#include "linalg.h"
#include "odeco.h"

struct odeco_euler_aligned {
	double theta;
	vec2 lambda;
	odeco_euler_aligned(double theta, const vec2& lambda);
	odeco_euler_aligned(const vec2& lambda);
};

struct odeco_mat_aligned {
	mat15 rot;
	vec2 lambda;

	odeco_mat_aligned(double theta, const vec2& lambda);
	odeco_mat_aligned(const vec2& lambda);
};

odeco_euler_aligned operator+(const odeco_euler_aligned& f0, const odeco_euler_aligned& f1);
odeco_euler_aligned operator*(double t, const odeco_euler_aligned& f);

vec15 ref_frame_coords_aligned(const vec3& lambda);
vec15 odeco_frame_coords_aligned(const odeco_euler_aligned& frame);
vec15 odeco_frame_coords_aligned(const odeco_mat_aligned& frame);

double loss(const vec15& y, const odeco_euler_aligned& frame);
double loss(const vec15& y, const odeco_mat_aligned& frame);

vec3 compute_gradient_aligned(const vec15& y, const odeco_euler_aligned& frame);
mat3 compute_hessian_aligned(const vec15& y, const odeco_euler_aligned& frame);

odeco_mat_aligned odeco_frame_project_aligned(const vec15& y, const vec3& d, double tol, int& num_iterations);
odeco_mat_aligned odeco_frame_project_aligned_newton(const vec15& y, double tol, int& num_iterations);