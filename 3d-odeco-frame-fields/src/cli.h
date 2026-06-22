#pragma once

#include <iostream>
#include "odeco.h"
#include "odeco_aligned.h"
#include "linalg.h"
#include <numbers>
#include <iomanip>
#include <chrono>
#include <sstream>

struct program_in {
	// 0: odeco_frame_project
	// 1: odeco_frame_project_aligned (to z axis)
	int mode;
	vec15 target;
	vec3 direction;
};