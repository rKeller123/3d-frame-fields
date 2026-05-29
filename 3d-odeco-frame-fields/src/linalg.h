#pragma once

#include <initializer_list>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <numbers>
#include <iostream>

// eigen type definitions

typedef Eigen::Matrix<double, 3, 1> vec3;
typedef Eigen::Matrix<double, 5, 1> vec5;
typedef Eigen::Matrix<double, 6, 1> vec6;
typedef Eigen::Matrix<double, 9, 1> vec9;
typedef Eigen::Matrix<double, 15, 1> vec15;

typedef Eigen::Matrix<double, 3, 3> mat3;
typedef Eigen::Matrix<double, 5, 5> mat5;
typedef Eigen::Matrix<double, 6, 6> mat6;
typedef Eigen::Matrix<double, 9, 9> mat9;
typedef Eigen::Matrix<double, 15, 15> mat15;
typedef Eigen::Matrix<double, 15, 3> mat15x3;
typedef Eigen::Matrix<double, 3, 15> mat3x15;

// rotations

mat5 rotation_z_5d(double gamma);
mat5 rotation_y_5d(double beta);
mat5 rotation_x_5d(double alpha);

mat9 rotation_z_9d(double gamma);
mat9 rotation_y_9d(double beta);
mat9 rotation_x_9d(double alpha);

mat15 rotation_15d(double alpha, double beta, double gamma);
mat15 rotation_15d(const vec3& theta);

mat3 rotate_rodrigues(const vec3& axis, double theta);
mat15 rotation_from_axis_angle(const vec3 & axis, double theta);

mat15 lie_x_15d();
mat15 lie_y_15d();
mat15 lie_z_15d();

mat15 rotation_15d_lie(double alpha, double beta, double gamma);

// vector operations

mat15 block_prod_mat(const mat15& a, const mat15& b);
vec15 block_prod_vec(const mat15& a, const vec15& b);

template<typename DerivedV0, typename DerivedV1>
inline double dot(const Eigen::MatrixBase<DerivedV0>& v0, const Eigen::MatrixBase<DerivedV1>& v1)
{
	return v0.dot(v1);
}

vec3 cross(const vec3& v0, const vec3& v1);

double norm(const vec15& v);

// matrix operations

mat15 exp(const mat15& m);

vec6 eivals(const mat6& m);
Eigen::LLT<mat6> llt(const mat6& m);

double clamped_log(double x);