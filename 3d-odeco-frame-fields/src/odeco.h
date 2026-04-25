#pragma once

#include "linalg.h"
#include <numbers>


vec15 reference_frame(double l_x, double l_y, double l_z);
vec15 odeco_frame(double alpha, double beta, double gamma, double l_x, double l_y, double l_z);

double objective_function(const vec15 y, double alpha, double beta, double gamma, double l_x, double l_y, double l_z);