// 3d-odeco-frame-fields.cpp : Defines the entry point for the application.

#include "3d-odeco-frame-fields.h"

using namespace std;

int parse(const string& input)
{
	try {
		return stoi(input);
	}
	catch (const invalid_argument&) {
		throw invalid_argument("parse: '" + input + "' is not a valid integer");
	}
	catch (const out_of_range&) {
		throw out_of_range("parse: '" + input + "' is out of int range");
	}
}

int main(int argc, char *argv[])
{
	int n;
	if (argc < 2) {
		cout << "No argument provided. Using default setting (6)" << endl << endl;
		n = 6;
	}
	else {
		n = parse(argv[1]);
	}

	double step = (numbers::pi + 1.08) / (2 * n);

	auto first_start = chrono::steady_clock::now();

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			for (int k = 0; k < n; k++) {
				odeco_mat target_frame(step * vec3(i, j, k), vec3(2.5, 7, 0.5));

				vec15 target_coords = odeco_frame_coords(target_frame);

				int num_iter;
				auto start = chrono::steady_clock::now();
				odeco_mat projected_frame = odeco_frame_project(target_coords, 0, 1e-9, 1e-5, num_iter);
				auto end = chrono::steady_clock::now();

				chrono::duration<double> time_elapsed = end - start;

				double proj_loss = loss(target_coords, projected_frame);

				cout << setprecision(4);
				cout << i << "i " << j << "j " << k << "k | " << time_elapsed.count() << "s | " << num_iter << "it | " << proj_loss << " loss" << endl;
			}
		}
	}

	auto last_end = chrono::steady_clock::now();

	chrono::duration<double> time_elapsed = last_end - first_start;

	cout << time_elapsed.count() << "s total" << endl << endl;

	return 0;
}
