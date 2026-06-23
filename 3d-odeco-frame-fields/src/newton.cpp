#include "newton.h"
#include "odeco.h"
#include "odeco_aligned.h"

double compute_step_size(const vec15& y, const odeco_euler& frame, const vec6& newton_step, double neg_newton_decr, double t_init, double alpha, double tau) {
    {
        odeco_euler dx_frame(newton_step.block<3, 1>(0, 0), newton_step.block<3, 1>(3, 0));

        double loss_y_frame = loss(y, frame);
        double alpha_dot = alpha * neg_newton_decr;

        double t = t_init;

        while (loss(y, frame + t * dx_frame) > loss_y_frame + t * alpha_dot) {
            t *= tau;
        }

        return t;
    }
}

double compute_step_size(const vec15& y, const odeco_euler_aligned& frame, const vec3& newton_step, double neg_newton_decr, double t_init, double alpha, double tau)
{
    odeco_euler_aligned dx_frame(newton_step(0), newton_step.template block<2, 1>(1, 0));

    double loss_y_frame = loss(y, frame);
    double alpha_dot = alpha * neg_newton_decr;

    double t = t_init;

    while (loss(y, frame + t * dx_frame) > loss_y_frame + t * alpha_dot) {
        t *= tau;
    }

    return t;
}

odeco_euler update_frame(const odeco_euler& frame, const vec6& newton_step, double t)
{
    odeco_euler newton_step_frame(newton_step.template block<3, 1>(0, 0), newton_step.template block<3, 1>(3, 0));
    odeco_euler update = frame + t * newton_step_frame;
    return update;
}

odeco_euler_aligned update_frame(const odeco_euler_aligned& frame, const vec3& newton_step, double t)
{
    odeco_euler_aligned newton_step_frame(newton_step(0), newton_step.template block<2, 1>(1, 0));
    odeco_euler_aligned update = frame + t * newton_step_frame;
    return update;
}