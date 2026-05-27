
#include "cli.h"

using namespace std;

program_in parse(const std::string& input)
{
	std::stringstream ss(input);
	vec15 values = vec15::Constant(0);
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

	return program_in{ .mode = mode, .target = values };
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
		cerr << "Usage: \"mode(0 or 1);y0;...;y15\"" << endl;
		return 1;
	}

	string input = argv[1];

	program_in in = parse(input);

	int mode = in.mode;
	vec15 target = in.target;

	switch (mode) {
	case 0:
	{
		int num_iter;
		odeco_mat proj = odeco_frame_project(target, 0, 1e-9, num_iter);
		cout << stringify(num_iter, odeco_frame_coords(proj)) << endl;
		break;
	}
	case 1:
	{
		vec15 proj_aligned = odeco_frame_project_aligned(target, vec3(0, 0, 1));
		cout << stringify(0, proj_aligned) << endl;
		break;
	}
	}

	return 0;
}