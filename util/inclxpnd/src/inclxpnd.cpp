#include <string>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

int parse(ifstream &input_file)
{	
	for (string line; getline (input_file, line); ) {
		stringstream tokenizer;
		tokenizer << line;
		
		string token;
		tokenizer >> token;
		
		if (token == "#include") {
			tokenizer >> token;
			if (token.length() < 2) return 2;

			auto ch = token.front();
			if (ch != '"' && ch != '<') return 3;

			ch = token.back();
			if (ch != '"' && ch != '>') return 3;

			string filename(token.begin() + 1, token.end() - 1);
			ifstream new_input_file(filename);
			
			if (new_input_file.good()) {
				parse(new_input_file);				
			} else {
				cout << "*** error: cannot include file: " << filename << endl;
			}
		} else {
			cout << line << endl;
		}
	}

	return 0;
}

int main(int argc, char* argv[])
{
	if (argc < 2) {
		cout << "Very simple file include expander. Outputs to stdout." << endl;
		cout << "Usage:" << endl;
		cout << "\tinclxpnd <file>" << endl;
		return 1;
	}

	ifstream input_file{ string(argv[1]) };
	if (!input_file.good()) return 1;

	return parse(input_file);
}