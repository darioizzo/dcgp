#define BOOST_TEST_MODULE dcgp_evaluation_perf
#include <boost/test/unit_test.hpp>
#include <boost/timer/timer.hpp>
#include <iostream>
#include <audi/gdual.hpp>

#include "../include/expression.hpp"
#include "../include/function_set.hpp"

using namespace audi;
using gdual_d = audi::gdual<double>;

void perform_evaluations(unsigned int in,
                  unsigned int out,
                  unsigned int rows,
                  unsigned int columns,
                  unsigned int levels_back,
                  unsigned int arity,
                  unsigned int N,
                  std::vector<dcgp::basis_function> function_set)
{
    // Random numbers engine
    std::default_random_engine re(123);
    // Instatiate the expression
    dcgp::expression ex(in, out, rows, columns, levels_back, arity, function_set, 123);
    // We create the input data upfront and we do not time it.
    std::vector<gdual_d> dumb(in);
    std::vector<std::vector<gdual_d> > in_num(N, dumb);

    for (auto j = 0u; j < N; ++j)
    {
        for (auto i = 0u; i < in; ++i)
        {
            auto value = std::uniform_real_distribution<double>(-1, 1)(re);
            in_num[j][i] = gdual_d(value, "x" + std::to_string(i), 1);
        }
    }

    std::cout << "Performing " << N << " evaluations, in:" << in << " out:" << out << " rows:" << rows << " columns:" << columns << std::endl;
    {
        boost::timer::auto_cpu_timer t;
        for (auto i = 0u; i < N; ++i)
        {
            ex(in_num[i]);
        }
    }
}

/// This torture test is passed whenever it completes. It is meant to check for
/// the code stability when large number of mutations are performed
BOOST_AUTO_TEST_CASE(evaluation_speed)
{
    unsigned int N = 1000;

    dcgp::function_set function_set1({"sum","diff","mul","div"});
    dcgp::stream(std::cout, "Function set ", function_set1(), "\n");
    perform_evaluations(2,4,2,3,4, 2, N, function_set1());
    perform_evaluations(2,4,10,10,11, 2, N, function_set1());
    perform_evaluations(2,4,20,20,21, 2, N, function_set1());
    perform_evaluations(1,1,1,100,101, 2, N, function_set1());
    perform_evaluations(1,1,2,100,101, 2, N, function_set1());
    perform_evaluations(1,1,3,100,101, 2, N, function_set1());

    dcgp::function_set function_set2({"sum","mul","sig"});
    dcgp::stream(std::cout, "Function set ", function_set2(), "\n");
    perform_evaluations(2,4,2,3,4, 2, N, function_set2());
    perform_evaluations(2,4,10,10,11, 2, N, function_set2());
    perform_evaluations(2,4,20,20,21, 2, N, function_set2());
    perform_evaluations(1,1,1,100,101, 2, N, function_set2());
    perform_evaluations(1,1,2,100,101, 2, N, function_set2());
    perform_evaluations(1,1,3,100,101, 2, N, function_set2());
}
