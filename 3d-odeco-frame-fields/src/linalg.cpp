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
    static const double sqrt_3 = sqrt(3.0);
    const double cos_beta = cos(beta);
    const double sin_beta = sin(beta);
    const double cos_2beta = cos(2.0 * beta);
    const double sin_2beta = sin(2.0 * beta);

    return mat5{
        { cos_beta, -sin_beta,  0, 0, 0},
        { sin_beta, cos_beta,  0, 0, 0},
        { 0, 0, 0.25 + 0.75 * cos_2beta, sqrt_3 * sin_2beta / 2, sqrt_3 * (1 - cos_2beta) / 4 },
        { 0, 0, -sqrt_3 * sin_2beta / 2, cos_2beta, sin_2beta / 2},
        { 0, 0, sqrt_3 * (1 - cos_2beta) / 4, -sin_2beta / 2, 0.75 + 0.25 * cos_2beta}
    };
}

mat5 rotation_x_5d(double alpha)
{
    static const double sqrt_3 = sqrt(3.0);
    const double cos_alpha = cos(alpha);
    const double sin_alpha = sin(alpha);
    const double cos_2alpha = cos(2.0 * alpha);
    const double sin_2alpha = sin(2.0 * alpha);

    return mat5{
        {cos_alpha, 0, 0, -sin_alpha, 0},
        {0, cos_2alpha, -sqrt_3 * sin_2alpha / 2, 0, -sin_2alpha / 2},
        {0, sqrt_3 * sin_2alpha / 2, 0.25 + 0.75 * cos_2alpha, 0, -sqrt_3 * (1 - cos_2alpha) / 4},
        {sin_alpha, 0, 0, cos_alpha, 0},
        {0, sin_2alpha / 2, -sqrt_3 * (1 - cos_2alpha) / 4, 0, 0.75 + 0.25 * cos_2alpha}
    };
}

mat9 rotation_z_9d(double gamma)
{
    const double cos_4gamma = cos(4.0 * gamma);
    const double cos_3gamma = cos(3.0 * gamma);
    const double cos_2gamma = cos(2.0 * gamma);
    const double cos_gamma = cos(gamma);

    const double sin_4gamma = sin(4.0 * gamma);
    const double sin_3gamma = sin(3.0 * gamma);
    const double sin_2gamma = sin(2.0 * gamma);
    const double sin_gamma = sin(gamma);

    return mat9{
        {cos_4gamma, 0, 0, 0, 0, 0, 0, 0, sin_4gamma},
        {0, cos_3gamma, 0, 0, 0, 0, 0, sin_3gamma, 0},
        {0, 0, cos_2gamma, 0, 0, 0, sin_2gamma, 0, 0},
        {0, 0, 0, cos_gamma, 0, sin_gamma, 0, 0, 0},
        {0, 0, 0, 0, 1, 0, 0, 0, 0},
        {0, 0, 0, -sin_gamma, 0, cos_gamma, 0, 0, 0},
        {0, 0, -sin_2gamma, 0, 0, 0, cos_2gamma, 0, 0},
        {0, -sin_3gamma, 0, 0, 0, 0, 0, cos_3gamma, 0},
        {-sin_4gamma, 0, 0, 0, 0, 0, 0, 0, cos_4gamma}
    };
}

mat9 rotation_y_9d(double beta)
{
    static const double sqrt_2 = sqrt(2.0);
    static const double sqrt_5 = sqrt(5.0);
    static const double sqrt_7 = sqrt(7.0);
    static const double sqrt_10 = sqrt(10.0);
    static const double sqrt_14 = sqrt(14.0);
    static const double sqrt_35 = sqrt(35.0);
    static const double sqrt_70 = sqrt(70.0);

    const double cos_beta = cos(beta);
    const double sin_beta = sin(beta);
    const double cos_2beta = cos(2.0 * beta);
    const double sin_2beta = sin(2.0 * beta);
    const double cos_3beta = cos(3.0 * beta);
    const double sin_3beta = sin(3.0 * beta);
    const double cos_4beta = cos(4.0 * beta);
    const double sin_4beta = sin(4.0 * beta);

    return mat9{
        { (1.0 / 2.0) * pow(cos_beta, 3) + (1.0 / 2.0) * cos_beta, -1.0 / 16.0 * sqrt_2 * (7 * sin_beta + 3 * sin_3beta), (1.0 / 8.0) * sqrt_7 * (cos_beta - cos_3beta), (1.0 / 16.0) * sqrt_14 * (-3 * sin_beta + sin_3beta), 0, 0, 0, 0, 0 },
        { (1.0 / 16.0) * sqrt_2 * (7 * sin_beta + 3 * sin_3beta), (7.0 / 16.0) * cos_beta + (9.0 / 16.0) * cos_3beta, (1.0 / 16.0) * sqrt_14 * (sin_beta - 3 * sin_3beta), (3.0 / 16.0) * sqrt_7 * (cos_beta - cos_3beta), 0, 0, 0, 0, 0 },
        { (1.0 / 8.0) * sqrt_7 * (cos_beta - cos_3beta), (1.0 / 16.0) * sqrt_14 * (-sin_beta + 3 * sin_3beta), (1.0 / 8.0) * cos_beta + (7.0 / 8.0) * cos_3beta, -1.0 / 16.0 * sqrt_2 * (3 * sin_beta + 7 * sin_3beta), 0, 0, 0, 0, 0 },
        { (1.0 / 16.0) * sqrt_14 * (3 * sin_beta - sin_3beta), (3.0 / 16.0) * sqrt_7 * (cos_beta - cos_3beta), (1.0 / 16.0) * sqrt_2 * (3 * sin_beta + 7 * sin_3beta), (9.0 / 16.0) * cos_beta + (7.0 / 16.0) * cos_3beta, 0, 0, 0, 0, 0 },
        { 0, 0, 0, 0, (35.0 / 8.0) * pow(sin_beta, 4) - 5 * pow(sin_beta, 2) + 1, (1.0 / 32.0) * sqrt_10 * (2 * sin_2beta + 7 * sin_4beta), (1.0 / 4.0) * sqrt_5 * (6 - 7 * pow(sin_beta, 2)) * pow(sin_beta, 2), (1.0 / 4.0) * sqrt_70 * pow(sin_beta, 3) * cos_beta, (1.0 / 8.0) * sqrt_35 * pow(sin_beta, 4) },
        { 0, 0, 0, 0, -1.0 / 32.0 * sqrt_10 * (2 * sin_2beta + 7 * sin_4beta), 7 * pow(sin_beta, 4) - 29.0 / 4.0 * pow(sin_beta, 2) + 1, (1.0 / 16.0) * sqrt_2 * (-2 * sin_2beta + 7 * sin_4beta), (1.0 / 8.0) * sqrt_7 * (cos_2beta - cos_4beta), (1.0 / 4.0) * sqrt_14 * pow(sin_beta, 3) * cos_beta },
        { 0, 0, 0, 0, (1.0 / 4.0) * sqrt_5 * (6 - 7 * pow(sin_beta, 2)) * pow(sin_beta, 2), (1.0 / 16.0) * sqrt_2 * (2 * sin_2beta - 7 * sin_4beta), (7.0 / 2.0) * pow(sin_beta, 4) - 4 * pow(sin_beta, 2) + 1, (1.0 / 2.0) * sqrt_14 * sin_beta * pow(cos_beta, 3), (1.0 / 4.0) * sqrt_7 * (1 - pow(cos_beta, 4)) },
        { 0, 0, 0, 0, -1.0 / 4.0 * sqrt_70 * pow(sin_beta, 3) * cos_beta, (1.0 / 8.0) * sqrt_7 * (cos_2beta - cos_4beta), -1.0 / 2.0 * sqrt_14 * sin_beta * pow(cos_beta, 3), pow(sin_beta, 4) - 11.0 / 4.0 * pow(sin_beta, 2) + 1, (1.0 / 32.0) * sqrt_2 * (14 * sin_2beta + sin_4beta) },
        { 0, 0, 0, 0, (1.0 / 8.0) * sqrt_35 * pow(sin_beta, 4), -1.0 / 4.0 * sqrt_14 * pow(sin_beta, 3) * cos_beta, (1.0 / 4.0) * sqrt_7 * (1 - pow(cos_beta, 4)), (1.0 / 32.0) * sqrt_2 * (-14 * sin_2beta - sin_4beta), (1.0 / 8.0) * pow(sin_beta, 4) - pow(sin_beta, 2) + 1 }
    };
}

mat9 rotation_x_9d(double alpha)
{
    static const double sqrt_2 = sqrt(2.0);
    static const double sqrt_5 = sqrt(5.0);
    static const double sqrt_7 = sqrt(7.0);
    static const double sqrt_10 = sqrt(10.0);
    static const double sqrt_14 = sqrt(14.0);
    static const double sqrt_35 = sqrt(35.0);
    static const double sqrt_70 = sqrt(70.0);

    const double cos_alpha = cos(alpha);
    const double sin_alpha = sin(alpha);
    const double cos_2alpha = cos(2.0 * alpha);
    const double sin_2alpha = sin(2.0 * alpha);
    const double cos_3alpha = cos(3.0 * alpha);
    const double sin_3alpha = sin(3.0 * alpha);
    const double cos_4alpha = cos(4.0 * alpha);
    const double sin_4alpha = sin(4.0 * alpha);


    return mat9{
        { (7.0 / 8.0) * cos_alpha + (1.0 / 8.0) * cos_3alpha, 0, -1.0 / 8.0 * sqrt_7 * cos_alpha + (1.0 / 8.0) * sqrt_7 * cos_3alpha, 0, 0, (3.0 / 16.0) * sqrt_14 * sin_alpha - 1.0 / 16.0 * sqrt_14 * sin_3alpha, 0, -7.0 / 16.0 * sqrt_2 * sin_alpha - 3.0 / 16.0 * sqrt_2 * sin_3alpha, 0 },
        { 0, (7.0 / 8.0) * cos_2alpha + (1.0 / 8.0) * cos_4alpha, 0, -1.0 / 8.0 * sqrt_7 * cos_2alpha + (1.0 / 8.0) * sqrt_7 * cos_4alpha, (1.0 / 16.0) * sqrt_70 * sin_2alpha - 1.0 / 32.0 * sqrt_70 * sin_4alpha, 0, -1.0 / 8.0 * sqrt_14 * sin_2alpha - 1.0 / 16.0 * sqrt_14 * sin_4alpha, 0, -7.0 / 16.0 * sqrt_2 * sin_2alpha - 1.0 / 32.0 * sqrt_2 * sin_4alpha },
        { -1.0 / 8.0 * sqrt_7 * cos_alpha + (1.0 / 8.0) * sqrt_7 * cos_3alpha, 0, (1.0 / 8.0) * cos_alpha + (7.0 / 8.0) * cos_3alpha, 0, 0, -3.0 / 16.0 * sqrt_2 * sin_alpha - 7.0 / 16.0 * sqrt_2 * sin_3alpha, 0, (1.0 / 16.0) * sqrt_14 * sin_alpha - 3.0 / 16.0 * sqrt_14 * sin_3alpha, 0 },
        { 0, -1.0 / 8.0 * sqrt_7 * cos_2alpha + (1.0 / 8.0) * sqrt_7 * cos_4alpha, 0, (1.0 / 8.0) * cos_2alpha + (7.0 / 8.0) * cos_4alpha, -1.0 / 16.0 * sqrt_10 * sin_2alpha - 7.0 / 32.0 * sqrt_10 * sin_4alpha, 0, (1.0 / 8.0) * sqrt_2 * sin_2alpha - 7.0 / 16.0 * sqrt_2 * sin_4alpha, 0, (1.0 / 16.0) * sqrt_14 * sin_2alpha - 1.0 / 32.0 * sqrt_14 * sin_4alpha },
        { 0, -1.0 / 16.0 * sqrt_70 * sin_2alpha + (1.0 / 32.0) * sqrt_70 * sin_4alpha, 0, (1.0 / 16.0) * sqrt_10 * sin_2alpha + (7.0 / 32.0) * sqrt_10 * sin_4alpha, (5.0 / 16.0) * cos_2alpha + (35.0 / 64.0) * cos_4alpha + 9.0 / 64.0, 0, -1.0 / 8.0 * sqrt_5 * cos_2alpha + (7.0 / 32.0) * sqrt_5 * cos_4alpha - 3.0 / 32.0 * sqrt_5, 0, -1.0 / 16.0 * sqrt_35 * cos_2alpha + (1.0 / 64.0) * sqrt_35 * cos_4alpha + (3.0 / 64.0) * sqrt_35 },
        { -3.0 / 16.0 * sqrt_14 * sin_alpha + (1.0 / 16.0) * sqrt_14 * sin_3alpha, 0, (3.0 / 16.0) * sqrt_2 * sin_alpha + (7.0 / 16.0) * sqrt_2 * sin_3alpha, 0, 0, (9.0 / 16.0) * cos_alpha + (7.0 / 16.0) * cos_3alpha, 0, -3.0 / 16.0 * sqrt_7 * cos_alpha + (3.0 / 16.0) * sqrt_7 * cos_3alpha, 0 },
        { 0, (1.0 / 8.0) * sqrt_14 * sin_2alpha + (1.0 / 16.0) * sqrt_14 * sin_4alpha, 0, -1.0 / 8.0 * sqrt_2 * sin_2alpha + (7.0 / 16.0) * sqrt_2 * sin_4alpha, -1.0 / 8.0 * sqrt_5 * cos_2alpha + (7.0 / 32.0) * sqrt_5 * cos_4alpha - 3.0 / 32.0 * sqrt_5, 0, (1.0 / 4.0) * cos_2alpha + (7.0 / 16.0) * cos_4alpha + 5.0 / 16.0, 0, (1.0 / 8.0) * sqrt_7 * cos_2alpha + (1.0 / 32.0) * sqrt_7 * cos_4alpha - 5.0 / 32.0 * sqrt_7 },
        { (7.0 / 16.0) * sqrt_2 * sin_alpha + (3.0 / 16.0) * sqrt_2 * sin_3alpha, 0, -1.0 / 16.0 * sqrt_14 * sin_alpha + (3.0 / 16.0) * sqrt_14 * sin_3alpha, 0, 0, -3.0 / 16.0 * sqrt_7 * cos_alpha + (3.0 / 16.0) * sqrt_7 * cos_3alpha, 0, (7.0 / 16.0) * cos_alpha + (9.0 / 16.0) * cos_3alpha, 0 },
        { 0, (7.0 / 16.0) * sqrt_2 * sin_2alpha + (1.0 / 32.0) * sqrt_2 * sin_4alpha, 0, -1.0 / 16.0 * sqrt_14 * sin_2alpha + (1.0 / 32.0) * sqrt_14 * sin_4alpha, -1.0 / 16.0 * sqrt_35 * cos_2alpha + (1.0 / 64.0) * sqrt_35 * cos_4alpha + (3.0 / 64.0) * sqrt_35, 0, (1.0 / 8.0) * sqrt_7 * cos_2alpha + (1.0 / 32.0) * sqrt_7 * cos_4alpha - 5.0 / 32.0 * sqrt_7, 0, (7.0 / 16.0) * cos_2alpha + (1.0 / 64.0) * cos_4alpha + 35.0 / 64.0 }
    };
}


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

mat3 rotate_rodrigues(const vec3& axis, double theta)
{
    mat3 nnT = axis * axis.transpose();
    mat3 cross = mat3{
        {0, -axis.z(), axis.y()},
        {axis.z(), 0, -axis.x()},
        {-axis.y(), axis.x(), 0}
    };
    return nnT + cos(theta) * (mat3::Identity() - nnT) + sin(theta) * cross;
}

mat15 rotation_from_axis_angle(const vec3& axis, double theta)
{
    mat3 rot3 = rotate_rodrigues(axis, theta);

    vec3 euler = rot3.canonicalEulerAngles(0, 1, 2);

    return rotation_15d(euler);
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

vec3 cross(const vec3& v0, const vec3& v1)
{
    return v0.cross(v1);
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

double clamped_log(double x)
{
    if (x <= 0) {
        return -INFINITY;
    }
    return log(x);
}
