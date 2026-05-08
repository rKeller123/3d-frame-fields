
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
	vec15 proj = odeco_frame_project(target, 0, 1e-9, num_iter);

	cout << stringify(num_iter, proj) << endl;

	return 0;
}