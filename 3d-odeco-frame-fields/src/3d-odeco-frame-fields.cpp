// 3d-odeco-frame-fields.cpp : Defines the entry point for the application.

#include "3d-odeco-frame-fields.h"

using namespace std;

int main()
{
	odeco_mat target_frame((numbers::pi / 2) * vec3(0.25, 0.75, 0.4), vec3(1.2, 2.4, 3.9));

	vec15 target_coords = odeco_frame_coords(target_frame);

	auto start = chrono::steady_clock::now();
	int num_iter;
	odeco_mat projected_frame = odeco_frame_project(target_coords, 5, 1e-6, num_iter);
	auto end = chrono::steady_clock::now();

	chrono::duration<double> time_elapsed = end - start;

	cout << "Projected in " << time_elapsed.count() << "s" << endl << endl;
	cout << "Projected in " << num_iter << " iterations" << endl << endl;

	double proj_loss = loss(target_coords, projected_frame);

	cout << "Proj - Target loss: " << proj_loss << endl << endl;

	int n = 6;

	double step = (numbers::pi + 1.08) / (2 * n);

	auto first_start = chrono::steady_clock::now();

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			for (int k = 0; k < n; k++) {
				target_frame = odeco_mat(step * vec3(i, j, k), vec3(2.5, 7, 0.5));

				target_coords = odeco_frame_coords(target_frame);

				start = chrono::steady_clock::now();
				projected_frame = odeco_frame_project(target_coords, 0, 1e-6, num_iter);
				end = chrono::steady_clock::now();

				time_elapsed = end - start;

				proj_loss = loss(target_coords, projected_frame);

				cout << setprecision(4);
				cout << i << "i " << j << "j " << k << "k | " << time_elapsed.count() << "s | " << num_iter << "it | " << proj_loss << " loss" << endl;
			}
		}
	}

	auto last_end = chrono::steady_clock::now();

	time_elapsed = last_end - first_start;

	cout << time_elapsed.count() << "s total" << endl << endl;

	return 0;
}
