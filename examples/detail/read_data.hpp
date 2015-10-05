// reading a text file
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include "../../src/std_overloads.hpp"

using namespace std;

// This reads data files in the csv format specified by Andrew James in his CGP-Library-V2.2
bool read_data (std::vector<std::vector<double> >& in, std::vector<std::vector<double> >& out, std::string filename) {
    string line;
    ifstream myfile (filename);
    if (myfile.is_open())
    {
        unsigned n,m,N;
        getline (myfile,line,',');
        {istringstream iss(line);
        vector<string> tokens{istream_iterator<string>{iss}, istream_iterator<string>{}};
        n = stoi(tokens[0]);}

        getline (myfile,line,',');
        {istringstream iss(line);
        vector<string> tokens{istream_iterator<string>{iss}, istream_iterator<string>{}};
        m = stoi(tokens[0]);}

        getline (myfile,line,',');
        {istringstream iss(line);
        vector<string> tokens{istream_iterator<string>{iss}, istream_iterator<string>{}};
        N = stoi(tokens[0]);}

        in.clear(); in.resize(N);
        out.clear(); out.resize(N);

        for (auto i = 0u; i < N; ++i) {
            for (auto j = 0u; j < n; ++j) {
                getline (myfile,line,',');
                istringstream iss(line);
                vector<string> tokens2{istream_iterator<string>{iss}, istream_iterator<string>{}};
                in[i].push_back(stod(tokens2[0]));
            }
            for (auto j = 0u; j < m; ++j) {
                getline (myfile,line,',');
                istringstream iss(line);
                vector<string> tokens2{istream_iterator<string>{iss}, istream_iterator<string>{}};
                out[i].push_back(stod(tokens2[0]));
            }
        }
        myfile.close();
        return true;
    }
    return false;
}

