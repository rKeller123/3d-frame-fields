#include "linalg.h"

// wigner matrices

const mat5 rot_x_pi_over_2_5d{
    {0, 0, 0, -1, 0},
    {0, -1, 0, 0, 0},
    {0, 0, -1.0 / 2, 0, -sqrt(3) / 2},
    {1, 0, 0, 0, 0},
    {0, 0, -sqrt(3) / 2, 0, 1.0 / 2}
};

const mat9 rot_x_pi_over_2_9d{
    {0, 0, 0, 0, 0, sqrt(14) / 4, 0, -sqrt(2) / 4, 0 },
    {0, -3.0 / 4, 0, sqrt(7) / 4, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, sqrt(2) / 4, 0, sqrt(14) / 4, 0},
    {0, sqrt(7) / 4, 0, 3.0 / 4, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 3.0 / 8, 0, sqrt(5) / 4, 0, sqrt(35) / 8},
    {-sqrt(14) / 4, 0, -sqrt(2) / 4, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, sqrt(5) / 4, 0, 1.0 / 2, 0, -sqrt(7) / 4},
    {sqrt(2) / 4, 0, -sqrt(14) / 4, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, sqrt(35) / 8, 0, -sqrt(7) / 4, 0, 1.0 / 8}
};

mat5 rotation_z_5d(double gamma)
{
    return mat5{
        {cos(2 * gamma), 0, 0, 0, sin(2 * gamma)},
        {0, cos(gamma), 0, sin(gamma), 0},
        {0, 0, 1, 0, 0},
        {0, -sin(gamma), 0, cos(gamma), 0},
        {-sin(2 * gamma), 0, 0, 0, cos(2 * gamma)}
    };
}

mat9 rotation_z_9d(double gamma)
{
    return mat9{
        {cos(4 * gamma), 0, 0, 0, 0, 0, 0, 0, sin(4 * gamma)},
        {0, cos(3 * gamma), 0, 0, 0, 0, 0 , sin(3 * gamma), 0},
        {0, 0, cos(2 * gamma), 0, 0, 0, sin(2 * gamma), 0, 0},
        {0, 0, 0, cos(gamma), 0, sin(gamma), 0, 0, 0},
        {0, 0, 0, 0, 1, 0, 0, 0, 0},
        {0, 0, 0, -sin(gamma), 0, cos(gamma), 0, 0, 0},
        {0, 0, -sin(2 * gamma), 0, 0, 0, cos(2 * gamma), 0, 0},
        {0, -sin(3 * gamma), 0, 0, 0, 0, 0, cos(3 * gamma), 0},
        {-sin(4 * gamma), 0, 0, 0, 0, 0, 0, 0, cos(4 * gamma)}
    };
}

// lie exponentiation and rotation

mat5 L_x_5d{
    {0, 0, 0, -1, 0},
    {0, 0, -sqrt(3), 0, -1},
    {0, sqrt(3), 0, 0, 0},
    {1, 0, 0, 0, 0},
    {0, 1, 0, 0, 0}
};

mat5 L_y_5d{
    {0, -1, 0, 0, 0},
    {1, 0, 0, 0, 0},
    {0, 0, 0, sqrt(3), 0},
    {0, 0, -sqrt(3), 0, 1},
    {0, 0, 0, -1, 0},
};

mat5 L_z_5d{
    {0, 0, 0, 0, 2},
    {0, 0, 0, 1, 0},
    {0, 0, 0, 0, 0},
    {0, -1, 0, 0, 0},
    {-2, 0, 0, 0, 0}
};

mat9 L_x_9d{
    {0, 0, 0, 0, 0, 0, 0, -sqrt(2), 0},
    {0, 0, 0, 0, 0, 0, -sqrt(7.0 / 2), 0, -sqrt(2)},
    {0, 0, 0, 0, 0, -3/sqrt(2), 0, -sqrt(7.0 / 2), 0},
    {0, 0, 0, 0, -sqrt(10), 0, -3 / sqrt(2), 0, 0},
    {0, 0, 0, sqrt(10), 0, 0, 0, 0, 0},
    {0, 0, 3/sqrt(2), 0, 0, 0, 0, 0, 0},
    {0, sqrt(7.0 / 2), 0, 3/sqrt(2), 0, 0, 0, 0, 0},
    {sqrt(2), 0, sqrt(7.0 / 2), 0, 0, 0, 0, 0, 0},
    {0, sqrt(2), 0, 0, 0, 0, 0, 0, 0},
};

mat9 L_y_9d{
    {0, -sqrt(2), 0, 0, 0, 0, 0, 0, 0},
    {sqrt(2), 0, -sqrt(7.0/2), 0, 0, 0, 0, 0, 0},
    {0, sqrt(7.0 / 2), 0, -3/sqrt(2), 0, 0, 0, 0, 0},
    {0, 0, 3/sqrt(2), 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, sqrt(10), 0, 0, 0},
    {0, 0, 0, 0, -sqrt(10), 0, 3 / sqrt(2), 0, 0},
    {0, 0, 0, 0, 0, -3/sqrt(2), 0, sqrt(7.0 / 2), 0},
    {0, 0, 0, 0, 0, 0, -sqrt(7.0 / 2), 0, sqrt(2)},
    {0, 0, 0, 0, 0, 0, 0, -sqrt(2), 0},
};

mat9 L_z_9d{
    {0, 0, 0, 0, 0, 0, 0, 0, 4},
    {0, 0, 0, 0, 0, 0, 0, 3, 0},
    {0, 0, 0, 0, 0, 0, 2, 0, 0},
    {0, 0, 0, 0, 0, 1, 0, 0, 0},
    {0, 0, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, -1, 0, 0, 0, 0, 0},
    {0, 0, -2, 0, 0, 0, 0, 0, 0},
    {0, -3, 0, 0, 0, 0, 0, 0, 0},
    {-4, 0, 0, 0, 0, 0, 0, 0, 0}
};

// methods

mat15 rotation_z_15d(double gamma)
{
    mat5 m5 = rotation_z_5d(gamma);
    mat9 m9 = rotation_z_9d(gamma);

    mat15 m = mat15::Zero();
    m(0, 0) = 1.0;
    m.block<5, 5>(1, 1) = m5;
    m.block<9, 9>(6, 6) = m9;

    return m;
}

mat15 rotation_y_15d(double beta)
{
    mat15 rot_x_pi_over_2_15d = mat15::Zero();
    rot_x_pi_over_2_15d(0, 0) = 1.0;
    rot_x_pi_over_2_15d.block<5, 5>(1, 1) = rot_x_pi_over_2_5d;
    rot_x_pi_over_2_15d.block<9, 9>(6, 6) = rot_x_pi_over_2_9d;

    return rot_x_pi_over_2_15d * rotation_z_15d(beta) * rot_x_pi_over_2_15d.transpose();
}

mat15 rotation_x_15d(double alpha)
{
    return rotation_y_15d(std::numbers::pi / 2).transpose() * rotation_z_15d(alpha) * rotation_y_15d(std::numbers::pi / 2);
}

mat15 rotation_15d(double alpha, double beta, double gamma)
{
    return rotation_x_15d(alpha) * rotation_y_15d(beta) * rotation_z_15d(gamma);
}

mat15 rotation_15d(const vec3& theta)
{
    return rotation_15d(theta.x(), theta.y(), theta.z());
}

mat15 lie_x_15d()
{
    mat15 m = mat15::Zero();

    m.block<5, 5>(1, 1) = L_x_5d;
    m.block<9, 9>(6, 6) = L_x_9d;

    return m;
}

mat15 lie_y_15d()
{
    mat15 m = mat15::Zero();

    m.block<5, 5>(1, 1) = L_y_5d;
    m.block<9, 9>(6, 6) = L_y_9d;

    return m;
}

mat15 lie_z_15d()
{
    mat15 m = mat15::Zero();

    m.block<5, 5>(1, 1) = L_z_5d;
    m.block<9, 9>(6, 6) = L_z_9d;

    return m;
}

mat15 rotation_15d_lie(double alpha, double beta, double gamma)
{
    return exp(alpha * lie_x_15d()) * exp(beta * lie_y_15d()) * exp(gamma * lie_z_15d());
}

double norm(const vec15& v)
{
    return v.norm();
}

mat15 exp(const mat15& m)
{
    return m.exp();
}

vec6 eivals(const mat6& m)
{
    Eigen::SelfAdjointEigenSolver<mat6> solver(m);

    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Eigenvalue decomposition failed");
    }
    return solver.eigenvalues();
}

Eigen::LLT<mat6> llt(const mat6& m)
{
    Eigen::LLT<mat6> llt_decomp(m);
    return llt_decomp;
}
