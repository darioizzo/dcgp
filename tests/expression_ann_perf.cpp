#include <random>
#define BOOST_TEST_MODULE dcgp_expression_ann_test
#include <algorithm>
#include <audi/back_compatibility.hpp>
#include <audi/io.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/timer/timer.hpp>

#include <dcgp/expression_ann.hpp>
#include <dcgp/kernel_set.hpp>

using namespace dcgp;

void perform_sgd(unsigned int rows, unsigned int columns, unsigned int levels_back, unsigned int arity, unsigned int N,
                 unsigned bs, std::vector<dcgp::kernel<double>> kernel_set)
{
    // Dimensions in and out are fixed
    unsigned in = 3u;
    unsigned out = 2u;
    // Random numbers
    std::default_random_engine rd(123);
    std::mt19937 gen{rd()};
    std::normal_distribution<> norm(0., 1.);

    // Instatiate the expression
    expression_ann<double> ex(in, out, rows, columns, levels_back, arity, kernel_set, 123);
    // We create the input data upfront and we do not time it.
    ex.randomise_weights();
    ex.randomise_biases();
    std::vector<double> in_dummy(in, 0.);
    std::vector<double> out_dummy(out, 0.);
    std::vector<std::vector<double>> data(N, in_dummy);
    std::vector<std::vector<double>> label(N, out_dummy);
    for (auto &item : data) {
        std::generate(item.begin(), item.end(), [&norm, &gen]() { return norm(gen); });
    }
    for (auto i = 0u; i < label.size(); ++i) {
        label[i][0] = 1. / 5. * std::cos(data[i][0] + data[i][1] + data[i][2]) - data[i][0] * data[i][1];
        label[i][1] = data[i][0] * data[i][1] * data[i][2];
    }

    std::cout << "One epoch of sgd:  rows:" << rows << " columns:" << columns << std::endl;
    {
        boost::timer::auto_cpu_timer t;
        ex.sgd(data, label, 0.01, bs, "MSE");
    }
}

/// This torture test is passed whenever it completes. It is meant to check for
/// the code stability when large number of mutations are performed
BOOST_AUTO_TEST_CASE(evaluation_speed)
{
    unsigned int N = 1024;

    dcgp::kernel_set<double> kernel_set1({"sig", "tanh", "ReLu"});
    dcgp::stream(std::cout, "Function set ", kernel_set1(), "\n");
    perform_sgd(100, 3, 1, 100, N, 32, kernel_set1());
    perform_sgd(100, 4, 1, 100, N, 32, kernel_set1());
    perform_sgd(100, 5, 1, 100, N, 32, kernel_set1());
    perform_sgd(100, 10, 1, 100, N, 32, kernel_set1());
    perform_sgd(100, 10, 1, 100, N, 256, kernel_set1());
}
