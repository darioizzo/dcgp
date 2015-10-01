
#define BOOST_TEST_MODULE dcgp_function_calls_perf
#include <boost/test/unit_test.hpp>
#include <boost/timer/timer.hpp>
#include <vector>

#include "../src/wrapped_functions.hpp"
#include "../src/basis_function.hpp" //my_fun_type

/*
double my_fun_call(const std::vector<double>& a, const std::vector<double>& b, unsigned int N)
{
    clock_t begin = clock();
    for (auto i = 0u; i < N; ++i)
    {
        dcgp::my_sum(a[i],b[i]);
    }

    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    return elapsed_secs;
}

double my_fun_call_indirect(const std::vector<double>& a, const std::vector<double>& b, unsigned int N)
{
    dcgp::function_set sum({"sum"});
    dcgp::expression ex(1, 1, 3, 3, 3, sum(), 123);
    clock_t begin = clock();
    for (auto i = 0u; i < N; ++i)
    {
        ex.get_f()[0].m_f(a[i],b[i]);
    }

    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    return elapsed_secs;
} */

/// We test the speed of evauating sum(a+b) calling my_fun_type and my_d_fun_type 
BOOST_AUTO_TEST_CASE(function_calls)
{
    // Number of evaluations to try
    unsigned int N = 10000000;
    // Random generators
    std::random_device rd;
    std::default_random_engine re(rd());

    // Generating the data set
    std::vector<double> a(N), b(N);
    std::vector<std::vector<double> > a_vector(N), b_vector(N);


    for (auto j = 0u; j < N; ++j)
    {
        a[j] = std::uniform_real_distribution<double>(-1, 1)(re);
        b[j] = std::uniform_real_distribution<double>(-1, 1)(re);

    }

    // Starting the test
    std::cout << "Testing " << N << " function calls to the sigmoid function" << std::endl;
    {
        boost::timer::auto_cpu_timer t; // Sets up a timer
        for (auto i = 0u; i < N; ++i)
        {
            dcgp::my_sig(a[i],b[i]);
        }
    }

    std::cout << "Testing " << N << " std::function calls to the sigmoid function" << std::endl;
    dcgp::my_fun_type my_sig2(dcgp::my_sig<double>);
    {
        boost::timer::auto_cpu_timer t; // Sets up a timer
        for (auto i = 0u; i < N; ++i)
        {
            my_sig2(a[i],b[i]);
        }
    }
}

