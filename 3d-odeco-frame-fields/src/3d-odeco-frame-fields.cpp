// 3d-odeco-frame-fields.cpp : Defines the entry point for the application.

#include "3d-odeco-frame-fields.h"

using namespace std;

struct bench_summary {
	int n;
	double mean;
	double variance;
	double min;
	double max;
	double sum;
	double sumsq;

	bench_summary() {
		n = 0;
		mean = 0.0;
		variance = 0.0;
		min = 0.0;
		max = 0.0;
		sum = 0.0;
		sumsq = 0.0;
	}

	void add(double x) {
		n++;
		sum += x;
		sumsq += x*x;

		mean = sum / n;
		variance = (sumsq / n) - (mean * mean);

		min = std::min(x, min);
		max = std::max(x, max);
	}

	std::string toString() const {
		return "mean=" + std::to_string(mean) +
			   ", std="  + std::to_string(sqrt(variance)) +
			   ", min="  + std::to_string(min) +
			   ", max="  + std::to_string(max);
	}
};

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

void compute_projections(int n)
{
	double step = (numbers::pi + 1.08) / (2 * n);

	auto first_start = chrono::steady_clock::now();

	bench_summary iter{};

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			for (int k = 0; k < n; k++) {
				odeco_mat target_frame(step * vec3(i, j, k), vec3(2.5, 7, 0.5));

				vec15 target_coords = odeco_frame_coords(target_frame);

				int num_iter;
				auto start = chrono::steady_clock::now();
				odeco_mat projected_frame = odeco_frame_project(target_coords, 0, 1e-9, num_iter);
				auto end = chrono::steady_clock::now();

				chrono::duration<double> time_elapsed = end - start;

				double proj_loss = loss(target_coords, projected_frame);

				cout << setprecision(4);
				cout << i << "i " << j << "j " << k << "k | " << time_elapsed.count() << "s | " << num_iter << "it | " << proj_loss << " loss" << endl;
				iter.add(num_iter);
			}
		}
	}

	auto last_end = chrono::steady_clock::now();

	chrono::duration<double> time_elapsed = last_end - first_start;

	cout << time_elapsed.count() << "s total" << endl << endl;
	cout << iter.toString() << endl;
}

void benchmark_seed_random(int n, int seed)
{
	mt19937 gen(seed);

	// can we normalize the target vector before projecting?
	uniform_real_distribution<double> dist(-5.0, 5.0);

	vec15 target;

	auto first_start = chrono::steady_clock::now();

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			for (int k = 0; k < n; k++) {

				for (int l = 0; l < 15; l++) {
					target(l) = dist(gen);
				}
				//target(0, 0) = 1.0;
				//target.block<5, 1>(1, 0).normalize();
				//target.block<9, 1>(6, 0).normalize();

				int num_iter_ref;
				auto start = chrono::steady_clock::now();
				odeco_mat projected_frame_ref = odeco_frame_project(target, 0, 1e-9, num_iter_ref);
				auto end = chrono::steady_clock::now();

				chrono::duration<double> time_elapsed_ref = end - start;

				double proj_loss_ref = loss(target, projected_frame_ref);

				int num_iter_seed;
				start = chrono::steady_clock::now();
				odeco_mat projected_frame_seed = odeco_frame_project(target, 5, 1e-9, num_iter_seed);
				end = chrono::steady_clock::now();

				chrono::duration<double> time_elapsed_seed = end - start;

				double proj_loss_seed = loss(target, projected_frame_seed);

				cout << setprecision(4);
				double time_speedup = time_elapsed_ref.count() / time_elapsed_seed.count();
				double iter_ratio = (double)num_iter_ref / num_iter_seed;
				double loss_ratio = proj_loss_ref / proj_loss_seed;

				// > 1 means using a grid searched seed is better
				cout << "Time speedup: " << time_speedup
					<< " | Iter ratio: " << iter_ratio
					<< " | Loss ratio: " << loss_ratio
					<< endl;

			}
		}
	}

	auto last_end = chrono::steady_clock::now();

	chrono::duration<double> time_elapsed = last_end - first_start;

	cout << time_elapsed.count() << "s total" << endl << endl;
}

void compute_aligned_projections(int n, const vec3& d)
{
	double step = (numbers::pi + 1.08) / (2 * n);

	bench_summary iter{};

	auto first_start = chrono::steady_clock::now();

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			for (int k = 0; k < n; k++) {
				odeco_mat target_frame(step * vec3(i, j, k), vec3(2.5, 7, 0.5));

				vec15 target_coords = odeco_frame_coords(target_frame);

				auto start = chrono::steady_clock::now();
				int num_iter;
				odeco_mat_aligned projected_frame = odeco_frame_project_aligned(target_coords, d, 1e-9, num_iter);
				auto end = chrono::steady_clock::now();

				chrono::duration<double> time_elapsed = end - start;

				double proj_loss = loss(target_coords, projected_frame);

				cout << setprecision(4);
				cout << i << "i " << j << "j " << k << "k | " << time_elapsed.count() << "s | " << num_iter << "it | " << proj_loss << " loss" << endl;
				iter.add(num_iter);
			}
		}
	}

	auto last_end = chrono::steady_clock::now();

	chrono::duration<double> time_elapsed = last_end - first_start;

	cout << time_elapsed.count() << "s total" << endl << endl;
	cout << iter.toString() << endl;
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

	compute_projections(n);
	// benchmark_seed_random(n, 42);
	// compute_aligned_projections(n, vec3(0, 1, 0));

	return 0;
}
