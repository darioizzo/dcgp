#include <iostream>
#include <iomanip>
#include <ctime>
#include <random>
#include "../src/dcgp.hpp"


double perform_evaluations(unsigned int in,
                  unsigned int out,
                  unsigned int rows,
                  unsigned int columns,
                  unsigned int levels_back,
                  unsigned int number_of_evaluations,
                  std::vector<dcgp::basis_function> function_set)
{
    // Random numbers engine
    std::default_random_engine re(123);
    // Instatiate the expression
    dcgp::expression ex(in, out, rows, columns, levels_back, function_set, 123);
    // We create the input data upfront and we do not time it.
    std::vector<double> dumb(in);
    std::vector<std::vector<double> > in_num(number_of_evaluations, dumb);

    for (auto j = 0u; j < number_of_evaluations; ++j)
    {
        for (auto i = 0u; i < in; ++i)
        {
            in_num[j][i] = std::uniform_real_distribution<double>(-1, 1)(re);
        }
    }

    clock_t begin = clock();
    for (auto i = 0u; i < number_of_evaluations; ++i)
    {
        //std::cout << dcgp::my_sum(in_num[i][0],in_num[i][1]) << "\n";
        //std::cout << ex.get_f()[0].m_f(in_num[i][0],in_num[i][1]) << std::endl;
        //std::cout << ex.get_f()[0].m_df({in_num[i][0]},{in_num[i][1]}) << std::endl;
        //dcgp::my_sum(in_num[i][0],in_num[i][1]);
        //ex.get_f()[0].m_f(in_num[i][0],in_num[i][1]);
        ex.get_f()[0].m_df({in_num[i][0]},{in_num[i][1]});
    }

    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    //std::cout << "In: " << in << ", Out: " << out << ", Rows: " << rows << ", Cols: " << columns << ", Levels-back: " << levels_back << ", Function set: " << ex.get_f() << std::endl;
    //std::cout << number_of_evaluations << " evaluations: " << elapsed_secs << " seconds" << std::endl;
    return elapsed_secs;
}

/// This torture test is passed whenever it completes. It is meant to check for
/// the code stability when large number of mutations are performed
int main() {
    double cum_time1=0;
    dcgp::function_set function_set1({"sum","diff","mul","div"});
    cum_time1+=perform_evaluations(2,4,2,3,4, 100000, function_set1());
    cum_time1+=perform_evaluations(2,4,10,10,11, 100000, function_set1());
    cum_time1+=perform_evaluations(2,4,20,20,21, 100000, function_set1());
    cum_time1+=perform_evaluations(1,1,1,100,101, 100000, function_set1());
    cum_time1+=perform_evaluations(1,1,2,100,101, 100000, function_set1());
    cum_time1+=perform_evaluations(1,1,3,100,101, 100000, function_set1());


    double cum_time2=0;
    dcgp::function_set function_set2({"sum","diff","sqrt","pow"});
    cum_time2+=perform_evaluations(2,4,2,3,4, 100000, function_set2());
    cum_time2+=perform_evaluations(2,4,10,10,11, 100000, function_set2());
    cum_time2+=perform_evaluations(2,4,20,20,21, 100000, function_set2());
    cum_time2+=perform_evaluations(1,1,1,100,101, 100000, function_set2());
    cum_time2+=perform_evaluations(1,1,2,100,101, 100000, function_set2());
    cum_time2+=perform_evaluations(1,1,3,100,101, 100000, function_set2());
    std::cout << "Cumulative time " << function_set1() << ": " << cum_time1 << " [s]" << std::endl;
    std::cout << "Cumulative time " << function_set2() << ": " << cum_time2 << " [s]" << std::endl;
    return 0;
}

