#define BOOST_TEST_MODULE dcgp_expression_ann_test
#include <boost/test/included/unit_test.hpp>

#include <algorithm>
#include <audi/back_compatibility.hpp>
#include <audi/io.hpp>
#include <boost/timer/timer.hpp>
#include <random>

#include <dcgp/expression_ann.hpp>
#include <dcgp/kernel_set.hpp>

using namespace dcgp;

void perform_sgd(unsigned int rows, unsigned int columns, unsigned int levels_back, const std::vector<unsigned> &arity,
                 unsigned int N, unsigned bs, std::vector<dcgp::kernel<double>> kernel_set, unsigned parallel)
{
    // Dimensions in and out are fixed
    unsigned in = 3u;
    unsigned out = 2u;
    // Random numbers
    std::default_random_engine rd(123);
    std::mt19937 gen{rd()};
    std::normal_distribution<> norm(0., 1.);

    // Instatiate the expression
    expression_ann ex(in, out, rows, columns, levels_back, arity, kernel_set, 123);
    // We create the input data upfront and we do not time it.
    ex.randomise_weights(0., 1., 123u);
    ex.randomise_biases(0., 1., 123u);
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
        ex.sgd(data, label, 0.01, bs, "MSE", parallel);
    }
}

/// This torture test is passed whenever it completes. It is meant to check for
/// the code stability when large number of mutations are performed
BOOST_AUTO_TEST_CASE(evaluation_speed)
{
    unsigned int N = 1024;
    dcgp::kernel_set<double> kernel_set1({"sig", "tanh", "ReLu", "ISRU", "ELU"});
    audi::print("Function set ", kernel_set1(), "\n");
    audi::print("Non parallel\n");
    perform_sgd(50, 3, 1, {3, 50, 20}, N, 32, kernel_set1(), false);
    perform_sgd(10, 10, 1, {100, 100, 100, 100, 100, 100, 100, 100, 100, 100}, N, 32u, kernel_set1(), false);
    perform_sgd(10, 100, 1, std::vector<unsigned>(100, 10), N, 32u, kernel_set1(), false);
    perform_sgd(100, 10, 1, {100, 10, 100, 10, 100, 10, 100, 10, 100, 10}, N, 32u, kernel_set1(), false);
    perform_sgd(100, 10, 1, {100, 100, 100, 100, 100, 100, 100, 100, 100, 100}, N, 32u, kernel_set1(), false);
    perform_sgd(200, 10, 1, {100, 100, 100, 100, 100, 100, 100, 100, 100, 100}, N, 32u, kernel_set1(), false);
    audi::print("Parallel\n");
    perform_sgd(50, 3, 1, {3, 50, 20}, N, 32u, kernel_set1(), 32u);
    perform_sgd(10, 10, 1, {100, 100, 100, 100, 100, 100, 100, 100, 100, 100}, N, 32u, kernel_set1(), 32u);
    perform_sgd(10, 100, 1, std::vector<unsigned>(100, 10), N, 32, kernel_set1(), 32u);
    perform_sgd(100, 10, 1, {100, 10, 100, 10, 100, 10, 100, 10, 100, 10}, N, 32u, kernel_set1(), 32u);
    perform_sgd(100, 10, 1, {100, 100, 100, 100, 100, 100, 100, 100, 100, 100}, N, 32u, kernel_set1(), 32u);
    perform_sgd(200, 10, 1, {100, 100, 100, 100, 100, 100, 100, 100, 100, 100}, N, 32u, kernel_set1(), 32u);
    audi::print("Increasing parallelism\n");
    perform_sgd(100, 10, 1, {100, 100, 100, 100, 100, 100, 100, 100, 100, 100}, N, 32u, kernel_set1(), false);
    perform_sgd(100, 10, 1, {100, 100, 100, 100, 100, 100, 100, 100, 100, 100}, N, 32u, kernel_set1(), 1u);
    perform_sgd(100, 10, 1, {100, 100, 100, 100, 100, 100, 100, 100, 100, 100}, N, 32u, kernel_set1(), 2u);
    perform_sgd(100, 10, 1, {100, 100, 100, 100, 100, 100, 100, 100, 100, 100}, N, 32u, kernel_set1(), 4u);
    perform_sgd(100, 10, 1, {100, 100, 100, 100, 100, 100, 100, 100, 100, 100}, N, 32u, kernel_set1(), 8u);
    perform_sgd(100, 10, 1, {100, 100, 100, 100, 100, 100, 100, 100, 100, 100}, N, 32u, kernel_set1(), 16u);
    perform_sgd(100, 10, 1, {100, 100, 100, 100, 100, 100, 100, 100, 100, 100}, N, 32u, kernel_set1(), 32u);
}
