#pragma once
#include "linalg.h"
#include "quaternion.h"

class quaternion;
inline quaternion operator*(const quaternion& lhs, const quaternion& rhs);

class quaternion {
public:
    quaternion(const double real, const vec3& im) {
        this->real_ = real;
        this->im_ = im;
    }

    static quaternion from_angle_axis(const double angle, const vec3& axis) {
        const double half_angle = angle / 2.0;
        return {
            std::cos(half_angle),
            std::sin(half_angle) * axis.normalized()
        };
    }

    static quaternion rotate_x(const double angle) {
        const double half_angle = angle / 2.0;
        return {std::cos(half_angle), {std::sin(half_angle), 0.0, 0.0}};
    }

    [[nodiscard]] static quaternion rotate_y(const double angle) noexcept
    {
        const auto half_angle = angle / 2.0;
        return {std::cos(half_angle), {0.0, std::sin(half_angle), 0.0}};
    }

    [[nodiscard]] static quaternion rotate_z(const double angle) noexcept
    {
        const auto half_angle = angle / 2.0;
        return {std::cos(half_angle), {0.0, 0.0, std::sin(half_angle)}};
    }

    static quaternion from_euler(const vec3& theta) {
        return rotate_x(theta.x()) * rotate_y(theta.y()) * rotate_z(theta.z());
    }

    static quaternion dir_to_z(const vec3& dir) {
        const vec3 unit_dir = dir.normalized();

        const vec3 axis = {unit_dir.y(), -unit_dir.x(), 0.0};
        const double s = std::sqrt(2.0 * (1.0 + unit_dir.z()));
        return {s / 2.0, axis / 2.0};

    }

    static quaternion z_to_dir(const vec3& dir) {
        return dir_to_z(dir).conjugate();
    }

    double real() const { return this->real_; }

    vec3 im() const { return this->im_; }

    quaternion conjugate() const { return {real(), -im()}; }

    vec3 to_euler() const {
        const auto qw = real();
        const auto qx = im().x();
        const auto qy = im().y();
        const auto qz = im().z();

        const auto asin_input = std::clamp(2.0 * (qx * qz + qw * qy), -1.0, 1.0);
        const auto y_angle = std::asin(asin_input);

        const auto [z_angle, x_angle] = [=]
        {
            if (std::abs(asin_input) >= 1.0 - 10.0 * std::numeric_limits<double>::epsilon()) [[unlikely]] {
                return std::make_pair(std::copysign(2.0, asin_input) * std::atan2(qx, qw),
                                      0.0);
            }
            return std::make_pair(std::atan2(-2.0 * (qx * qy - qz * qw),
                                             qw * qw + qx * qx - qy * qy - qz * qz),
                                  std::atan2(-2.0 * (qy * qz - qx * qw),
                                             qw * qw - qx * qx - qy * qy + qz * qz));
        }();

        return {x_angle, y_angle, z_angle};
    }

    std::pair<double, vec3> to_angle_axis() const {
        const auto qw = real();
        const auto qx = im().x();
        const auto qy = im().y();
        const auto qz = im().z();
        
        const auto angle = 2 * std::acos(qw);
        const auto x = qx / std::sqrt(1 - qw * qw);
        const auto y = qy / std::sqrt(1 - qw * qw);
        const auto z = qz / std::sqrt(1 - qw * qw);

        return std::make_pair(angle, vec3(x, y, z));
    }

private:
    double real_;
    vec3 im_{};

};

inline quaternion operator*(const quaternion& lhs, const quaternion& rhs)
{
    return {lhs.real() * rhs.real() - lhs.im().dot(rhs.im()),
            lhs.real() * rhs.im() + rhs.real() * lhs.im() + lhs.im().cross(rhs.im())};
}
