// reading a text file
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

// This reads data files in the csv format specified by Andrew James in his CGP-Library-V2.2
bool read_data (std::vector<std::vector<double> >& in, std::vector<std::vector<double> >& out, std::string filename) {
    std::string line;
    std::ifstream myfile (filename);
    if (myfile.is_open())
    {
        unsigned n,m,N;
        getline (myfile,line,',');
        {
        std::istringstream iss(line);
        std::vector<std::string> tokens{std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>{}};
        n = stoi(tokens[0]);}

        getline (myfile,line,',');
        {
        std::istringstream iss(line);
        std::vector<std::string> tokens{std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>{}};
        m = std::stoi(tokens[0]);}

        getline (myfile,line,',');
        {
        std::istringstream iss(line);
        std::vector<std::string> tokens{std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>{}};
        N = std::stoi(tokens[0]);}

        in.clear(); in.resize(N);
        out.clear(); out.resize(N);

        for (auto i = 0u; i < N; ++i) {
            for (auto j = 0u; j < n; ++j) {
                getline (myfile,line,',');
                std::istringstream iss(line);
                std::vector<std::string> tokens2{std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>{}};
                in[i].push_back(std::stod(tokens2[0]));
            }
            for (auto j = 0u; j < m; ++j) {
                getline (myfile,line,',');
                std::istringstream iss(line);
                std::vector<std::string> tokens2{std::istream_iterator<std::string>{iss}, std::istream_iterator<std::string>{}};
                out[i].push_back(std::stod(tokens2[0]));
            }
        }
        myfile.close();
        return true;
    }
    return false;
}
