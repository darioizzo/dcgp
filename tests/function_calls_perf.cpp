
#define BOOST_TEST_MODULE dcgp_function_calls_perf
#include <boost/test/unit_test.hpp>
#include <boost/timer/timer.hpp>
#include <vector>

#include "../include/wrapped_functions.hpp"
#include "../include/basis_function.hpp" //my_fun_type
#include "../include/function_set.hpp"
#include "../include/expression.hpp"


// We test the speed of evauating sig(a+b) calling
// the function directly, via an std::function or a minimal d-CGP expression

BOOST_AUTO_TEST_CASE(function_calls)
{
    // Number of evaluations to try
    unsigned int N = 100000;
    // Random generators
    std::random_device rd;
    std::default_random_engine re(rd());

    // Generating the data set
    std::vector<double> a(N), b(N);
    std::vector<std::vector<double> > ab_vector(N);

    for (auto j = 0u; j < N; ++j)
    {
        a[j] = std::uniform_real_distribution<double>(-1, 1)(re);
        b[j] = std::uniform_real_distribution<double>(-1, 1)(re);
        ab_vector[j] = {a[j], b[j]};
    }

    // Starting the test
    std::cout << "Testing " << N << " function calls to the sigmoid function" << std::endl;
    {
        boost::timer::auto_cpu_timer t; // Sets up a timer
        for (auto i = 0u; i < N; ++i)
        {
            dcgp::my_sig<double>({a[i],b[i]});
        }
    }

    std::cout << "Testing " << N << " std::function calls to the sigmoid function" << std::endl;
    dcgp::my_fun_type my_sig2(dcgp::my_sig<double>);
    {
        boost::timer::auto_cpu_timer t; // Sets up a timer
        for (auto i = 0u; i < N; ++i)
        {
            my_sig2({a[i],b[i]});
        }
    }

    std::cout << "Testing " << N << " std::function calls to the sigmoid function via dcgp::expression" << std::endl;
    dcgp::function_set sigmoid_set({"sig"});
    dcgp::expression ex(2,1,1,1,1,2,sigmoid_set(),0);
    ex.set({0,0,1,2});
    {
        boost::timer::auto_cpu_timer t; // Sets up a timer
        for (auto i = 0u; i < N; ++i)
        {
            ex(ab_vector[i]);
        }
    }
}
