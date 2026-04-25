#include <initializer_list>
#include <Eigen/Dense>
#include <numbers>

// eigen type definitions

typedef Eigen::Matrix<double, 3, 1> vec3;
typedef Eigen::Matrix<double, 5, 1> vec5;
typedef Eigen::Matrix<double, 9, 1> vec9;
typedef Eigen::Matrix<double, 15, 1> vec15;

typedef Eigen::Matrix<double, 3, 3> mat3;
typedef Eigen::Matrix<double, 5, 5> mat5;
typedef Eigen::Matrix<double, 9, 9> mat9;
typedef Eigen::Matrix<double, 15, 15> mat15;


// rotations

mat5 rotation_z_5d(double gamma);

mat9 rotation_z_9d(double gamma);

mat15 rotation_z_15d(double gamma);
mat15 rotation_y_15d(double beta);
mat15 rotation_x_15d(double alpha);
mat15 rotation_15d(double alpha, double beta, double gamma);

// vector operations


double dot(vec15 v0, vec15 v1);
double norm(vec15 v);