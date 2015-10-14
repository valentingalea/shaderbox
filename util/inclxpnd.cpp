#include <string>
#include <fstream>
#include <iostream>
#include <sstream>

int main ()
{
	using namespace std;
	
	ifstream input_file("../src/app_2d.h");
	if (input_file.bad ()) return 1;
	
	for (string line; getline (input_file, line); ) {
		stringstream tokenizer;
		tokenizer << line;
		
		string token;
		tokenizer >> token;
		
		if (token == "#include") {
			cout << "++++++++++" << endl;
		} else {
			cout << line << endl;
		}
	}
	
	return 0;
}