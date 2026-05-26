
#include "cli.h"

using namespace std;

vec15 parse(const string& input)
{
	stringstream ss(input);
	vec15 values = vec15::Constant(0);
	string item;

	for (int i = 0; i < 15; i++) {
		if (!getline(ss, item, ';')) {
			throw runtime_error("Not enough values in input string");
		}
		try {
			values(i) = std::stod(item);
		}
		catch (...) {
			throw std::runtime_error("Invalid number: " + item);
		}
		
	}
	return values;
}

string stringify(int num_iter, const vec15& v)
{
	string result = to_string(num_iter) + ";";
	
	for (int i = 0; i < 15; i++) {
		result.append(to_string(v(i)));

		if (i < 14) {
			result.append(";");
		}
	}
	
	return result;
}

int main(int argc, char *argv[])
{
	if (argc < 2) {
		cerr << "Usage: Provide a string of 15 doubles seperated by ';'" << endl;
		return 1;
	}

	string input = argv[1];

	vec15 target = parse(input);

	int num_iter;
	odeco_mat proj = odeco_frame_project(target, 0, 1e-9, num_iter);

	static const double M_PI = numbers::pi;

	vec15 z_aligned = vec15(
		(7.0 / 10.0) * target(0) + (1.0 / 10.0) * target(10) - 1.0 / 5.0 * sqrt(5) * target(3) + (1.0 / 3.0) * sqrt(M_PI),
		(27.0 / 28.0) * target(1) - 3.0 / 28.0 * sqrt(3) * target(8),
		0,
		(1.0 / 105.0) * sqrt(5) * (-21 * target(0) - 3 * target(10) + 6 * sqrt(5) * target(3) + 26 * sqrt(M_PI)),
		0,
		-3.0 / 28.0 * sqrt(3) * target(12) + (27.0 / 28.0) * target(5),
		target(6),
		0,
		-3.0 / 28.0 * sqrt(3) * target(1) + (1.0 / 28.0) * target(8),
		0,
		(1.0 / 10.0) * target(0) + (1.0 / 70.0) * target(10) - 1.0 / 35.0 * sqrt(5) * target(3) + (1.0 / 7.0) * sqrt(M_PI),
		0,
		(1.0 / 28.0) * target(12) - 3.0 / 28.0 * sqrt(3) * target(5),
		0,
		target(14)
	);

	cout << stringify(num_iter, z_aligned) << endl;

	return 0;
}