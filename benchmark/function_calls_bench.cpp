#include <boost/timer/timer.hpp>
#include <vector>

#include <dcgp/expression.hpp>
#include <dcgp/function.hpp>
#include <dcgp/kernel_set.hpp>
#include <dcgp/wrapped_functions.hpp>

// We test the speed of evauating sig(a+b) calling
// the function directly, via a function or a minimal d-CGP expression

int main()
{
    // Number of evaluations to try
    unsigned int N = 100000;
    // Random generators
    std::random_device rd;
    std::default_random_engine re(rd());

    // Generating the data set
    std::vector<double> a(N), b(N);
    std::vector<std::vector<double>> ab_vector(N);

    for (auto j = 0u; j < N; ++j) {
        a[j] = std::uniform_real_distribution<double>(-1, 1)(re);
        b[j] = std::uniform_real_distribution<double>(-1, 1)(re);
        ab_vector[j] = {a[j], b[j]};
    }

    // Starting the test
    std::cout << "Testing " << N << " function calls to the sigmoid function" << std::endl;
    {
        boost::timer::auto_cpu_timer t; // Sets up a timer
        for (auto i = 0u; i < N; ++i) {
            dcgp::my_sig<double>({a[i], b[i]});
        }
    }

    std::cout << "Testing " << N << " function calls to the sigmoid function" << std::endl;
    dcgp::function<double(const std::vector<double> &)> my_sig2(dcgp::my_sig<double>);
    {
        boost::timer::auto_cpu_timer t; // Sets up a timer
        for (auto i = 0u; i < N; ++i) {
            my_sig2({a[i], b[i]});
        }
    }

    std::cout << "Testing " << N << " function calls to the sigmoid function via dcgp::expression" << std::endl;
    dcgp::kernel_set<double> only_one_sigmoid({"sig"});
    dcgp::expression<double> ex(2, 1, 1, 1, 1, 2, only_one_sigmoid(), 0u, 0u);
    ex.set({0, 0, 1, 2});
    {
        boost::timer::auto_cpu_timer t; // Sets up a timer
        for (auto i = 0u; i < N; ++i) {
            ex(ab_vector[i]);
        }
    }
    return 0;
}
