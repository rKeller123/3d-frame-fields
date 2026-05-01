// 3d-odeco-frame-fields.cpp : Defines the entry point for the application.

#include "3d-odeco-frame-fields.h"

using namespace std;

int main()
{
	odeco_frame_description target_frame = {
		.theta = vec3(numbers::pi / 2, 3 * numbers::pi / 4, 3 * numbers::pi / 24),
		.lambda = vec3(3, 6, 9)
	};

	vec15 target = odeco_frame(target_frame);

	odeco_frame_description projected_description = odeco_frame_project(target, 5, 1e-6);

	cout << "Target frame description:    " << endl << target_frame.theta << endl << endl << target_frame.lambda << endl << endl;
	cout << "Projected frame description: " << endl << projected_description.theta << endl << endl << projected_description.lambda << endl << endl;

	double proj_loss = loss(target, projected_description);

	cout << "Proj - Target Loss: " << proj_loss << endl << endl;

	return 0;
}
