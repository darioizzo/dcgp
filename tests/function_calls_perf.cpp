#include <iostream>
#include <iomanip>
#include <ctime>
#include <random>
#include "../src/dcgp.h"


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

double d_my_fun_call(const std::vector<std::vector<double> >& a, const std::vector<std::vector<double> >& b, unsigned int N)
{
    clock_t begin = clock();
    for (auto i = 0u; i < N; ++i)
    {
        dcgp::d_my_sum(a[i],b[i]);
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
}

double d_my_fun_call_indirect(const std::vector<std::vector<double> >& a, const std::vector<std::vector<double> >& b, unsigned int N)
{
    dcgp::function_set sum({"sum"});
    dcgp::expression ex(1, 1, 3, 3, 3, sum(), 123);
    clock_t begin = clock();
    for (auto i = 0u; i < N; ++i)
    {
        ex.get_f()[0].m_df(a[i],b[i]);
    }

    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    return elapsed_secs;
}

/// We test the speed of evauating sum(a+b) calling my_fun_type and my_d_fun_type 
int main() {
    // Number of evaluations tried
    unsigned int N = 1000000;
    // Generating the data set
    std::vector<double> a(N), b(N);
    std::vector<std::vector<double> > a_vector(N), b_vector(N);
    std::default_random_engine re(123);

    for (auto j = 0u; j < N; ++j)
    {
        a[j] = std::uniform_real_distribution<double>(-1, 1)(re);
        b[j] = std::uniform_real_distribution<double>(-1, 1)(re);
        a_vector[j] = {a[j]};
        b_vector[j] = {b[j]};
    }

    double time_my_fun=my_fun_call(a, b, N);
    double time_my_fun_indirect=my_fun_call_indirect(a, b, N);
    double time_d_my_fun=d_my_fun_call(a_vector, b_vector, N);
    double time_d_my_fun_indirect=d_my_fun_call_indirect(a_vector, b_vector, N);
    std::cout << "Cumulative time my_fun: " << time_my_fun << " [s]" << std::endl;
    std::cout << "Cumulative time my_fun_indirect: " << time_my_fun_indirect << " [s]" << std::endl;
    std::cout << "Cumulative time d_my_fun: " << time_d_my_fun << " [s]" << std::endl;
    std::cout << "Cumulative time d_my_fun_indirect: " << time_d_my_fun_indirect << " [s]" << std::endl;
    return 0;
}

