#pragma once

#include "linalg.h"
#include "odeco.h"

struct odeco_euler_aligned {
	double theta;
	vec2 lambda;
	odeco_euler_aligned(double theta, const vec2& lambda);
	odeco_euler_aligned(const vec2& lambda);
};

struct odeco_quat_aligned {
	quaternion rot;
	vec2 lambda;

	odeco_quat_aligned(const quaternion& rot, const vec2& lambda);
	odeco_quat_aligned(double theta, const vec2& lambda);
	odeco_quat_aligned(const vec2& lambda);
};

odeco_euler_aligned operator+(const odeco_euler_aligned& f0, const odeco_euler_aligned& f1);
odeco_euler_aligned operator*(double t, const odeco_euler_aligned& f);

vec15 ref_frame_coords_aligned(const vec2& lambda);
vec15 odeco_frame_coords(const odeco_euler_aligned& frame);
vec15 odeco_frame_coords(const odeco_quat_aligned& frame);

double loss(const vec15& y, const odeco_euler_aligned& frame);
double loss(const vec15& y, const odeco_quat_aligned& frame);

vec3 compute_gradient(const vec15& y, const odeco_euler_aligned& frame);
mat3 compute_hessian(const vec15& y, const odeco_euler_aligned& frame);

odeco_quat_aligned odeco_frame_project_aligned(const vec15& y, const vec3& d, double tol, int& num_iterations);
odeco_quat_aligned odeco_frame_project_aligned_newton(const vec15& y, double tol, int& num_iterations);