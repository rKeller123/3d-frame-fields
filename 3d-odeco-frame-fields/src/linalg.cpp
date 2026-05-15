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
    {0, 0, 0, 0, 0, -3 / sqrt(2), 0, -sqrt(7.0 / 2), 0},
    {0, 0, 0, 0, -sqrt(10), 0, -3 / sqrt(2), 0, 0},
    {0, 0, 0, sqrt(10), 0, 0, 0, 0, 0},
    {0, 0, 3 / sqrt(2), 0, 0, 0, 0, 0, 0},
    {0, sqrt(7.0 / 2), 0, 3 / sqrt(2), 0, 0, 0, 0, 0},
    {sqrt(2), 0, sqrt(7.0 / 2), 0, 0, 0, 0, 0, 0},
    {0, sqrt(2), 0, 0, 0, 0, 0, 0, 0},
};

mat9 L_y_9d{
    {0, -sqrt(2), 0, 0, 0, 0, 0, 0, 0},
    {sqrt(2), 0, -sqrt(7.0 / 2), 0, 0, 0, 0, 0, 0},
    {0, sqrt(7.0 / 2), 0, -3 / sqrt(2), 0, 0, 0, 0, 0},
    {0, 0, 3 / sqrt(2), 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, sqrt(10), 0, 0, 0},
    {0, 0, 0, 0, -sqrt(10), 0, 3 / sqrt(2), 0, 0},
    {0, 0, 0, 0, 0, -3 / sqrt(2), 0, sqrt(7.0 / 2), 0},
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

mat5 rotation_y_5d(double beta)
{
    static const mat5 rot_x_pi_over_2_5d_T = rot_x_pi_over_2_5d.transpose();
    mat5 tmp;
    tmp.noalias() = rotation_z_5d(beta) * rot_x_pi_over_2_5d_T;
    return rot_x_pi_over_2_5d * tmp;
}

mat5 rotation_x_5d(double alpha)
{
    static const mat5 rot_y_pi_over_2_5d = rotation_y_5d(std::numbers::pi / 2);
    static const mat5 rot_y_pi_over_2_5d_T = rot_y_pi_over_2_5d.transpose();
    mat5 tmp;
    tmp.noalias() = rotation_z_5d(alpha) * rot_y_pi_over_2_5d;
    return rot_y_pi_over_2_5d_T * tmp;
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

mat9 rotation_y_9d(double beta)
{
    static const mat9 rot_x_pi_over_2_9d_T = rot_x_pi_over_2_9d.transpose();
    mat9 tmp;
    tmp.noalias() = rotation_z_9d(beta) * rot_x_pi_over_2_9d_T;
    return rot_x_pi_over_2_9d * tmp;
}

mat9 rotation_x_9d(double alpha)
{
    static const mat9 rot_y_pi_over_2_9d = rotation_y_9d(std::numbers::pi / 2);
    static const mat9 rot_y_pi_over_2_9d_T = rot_y_pi_over_2_9d.transpose();
    mat9 tmp;
    tmp.noalias() = rotation_z_9d(alpha) * rot_y_pi_over_2_9d;
    return rot_y_pi_over_2_9d_T * tmp;
}

// methods



mat15 rotation_15d(double alpha, double beta, double gamma)
{
    mat5 rot_x_5d = rotation_x_5d(alpha);
    mat5 rot_y_5d = rotation_y_5d(beta);
    mat5 rot_z_5d = rotation_z_5d(gamma);

    mat9 rot_x_9d = rotation_x_9d(alpha);
    mat9 rot_y_9d = rotation_y_9d(beta);
    mat9 rot_z_9d = rotation_z_9d(gamma);

    mat15 rot = mat15::Zero();

    // leverage block matrix structure for more efficient product
    rot(0, 0) = 1.0;

    mat5 tmp5;
    tmp5.noalias() = rot_y_5d * rot_z_5d;
    mat5 res5;
    res5.noalias() = rot_x_5d * tmp5;
    rot.block<5, 5>(1, 1) = res5;

    mat9 tmp9;
    tmp9.noalias() = rot_y_9d * rot_z_9d;
    mat9 res9;
    res9.noalias() = rot_x_9d * tmp9;
    rot.block<9, 9>(6, 6) = res9;
  
    return rot;
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

mat15 block_prod_mat(const mat15& a, const mat15& b)
{
    mat15 res = mat15::Zero();

    res(0, 0) = a(0, 0) * b(0, 0);
    res.block<5, 5>(1, 1) = a.block<5, 5>(1, 1) * b.block<5, 5>(1, 1);
    res.block<9, 9>(6, 6) = a.block<9, 9>(6, 6) * b.block<9, 9>(6, 6);

    return res;
}

vec15 block_prod_vec(const mat15& a, const vec15& b)
{
    vec15 res = vec15::Zero();

    res(0) = a(0, 0) * b(0);
    res.block<5, 1>(1, 0) = a.block<5, 5>(1, 1) * b.block<5, 1>(1, 0);
    res.block<9, 1>(6, 0) = a.block<9, 9>(6, 6) * b.block<9, 1>(6, 0);
    
    return res;
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
