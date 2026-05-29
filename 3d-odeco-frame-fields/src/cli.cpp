
#include "cli.h"

using namespace std;

program_in parse(const std::string& input)
{
	std::stringstream ss(input);
	vec15 values = vec15::Zero();
	vec3 direction = vec3::Zero();
	std::string item;

	int mode;

	if (!std::getline(ss, item, ';')) {
		std::cerr << "Error: Not enough values in input string\n";
		std::exit(1);
	}

	try {
		mode = std::stoi(item);
	}
	catch (const std::exception&) {
		std::cerr << "Error: Invalid mode: " << item << "\n";
		std::exit(1);
	}

	if (mode != 0 && mode != 1) {
		std::cerr << "Error: Invalid mode: " << mode << "\n";
		std::exit(1);
	}

	if (mode == 1) {
		for (int i = 0; i < 3; i++) {
			if (!std::getline(ss, item, ';')) {
				std::cerr << "Error: Not enough values for align direction in input string\n";
				std::exit(1);
			}
			try {
				direction(i) = std::stod(item);
			}
			catch (const std::exception&) {
				std::cerr << "Error: Invalid number: " << item << "\n";
				std::exit(1);
			}
		}
	}

	for (int i = 0; i < 15; i++) {
		if (!std::getline(ss, item, ';')) {
			std::cerr << "Error: Not enough values in input string\n";
			std::exit(1);
		}

		try {
			values(i) = std::stod(item);
		}
		catch (const std::exception&) {
			std::cerr << "Error: Invalid number: " << item << "\n";
			std::exit(1);
		}
	}

	direction.normalize();
	return program_in{ .mode = mode, .target = values, .direction = direction };
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
		cerr << "Usage: \"mode(0 or 1);(d0;d1;d2 if mode = 1);y0;...;y15\"" << endl;
		return 1;
	}

	string input = argv[1];

	program_in in = parse(input);

	int mode = in.mode;
	vec15 target = in.target;
	vec3 direction = in.direction;

	int num_iter;
	switch (mode) {
	case 0:
	{
		odeco_mat proj = odeco_frame_project(target, 0, 1e-9, num_iter);
		cout << stringify(num_iter, odeco_frame_coords(proj)) << endl;
		break;
	}
	case 1:
	{
		odeco_mat proj_aligned = odeco_frame_project_aligned(target, direction, 0, 1e-9, num_iter);
		cout << stringify(num_iter, odeco_frame_coords(proj_aligned)) << endl;
		break;
	}
	}

	return 0;
}