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

void perform_evaluations(unsigned int in, unsigned int out, unsigned int rows, unsigned int columns,
                         unsigned int levels_back, unsigned int arity, unsigned int N,
                         std::vector<dcgp::kernel<double>> kernel_set)
{
    // Random numbers engine
    std::default_random_engine rd(123);
    std::mt19937 gen{rd()};
    std::normal_distribution<> norm(0., 1.);

    // Instatiate the expression
    expression_ann<double> ex(in, out, rows, columns, levels_back, arity, kernel_set, 123);
    // We create the input data upfront and we do not time it.
    ex.randomise_weights();
    ex.randomise_biases();
    std::vector<std::vector<double>> data(N, std::vector<double>(in, 0.));
    std::vector<std::vector<double>> label(N, std::vector<double>(out, 0.));
    for (auto &item : data) {
        std::generate(item.begin(), item.end(), [&norm, &gen]() { return norm(gen); });
    }
    for (auto i = 0u; i < label.size(); ++i) {
        label[i][0] = 1. / 5. * std::cos(data[i][0] + data[i][1] + data[i][2]) - data[i][0] * data[i][1];
        label[i][1] = data[i][0] * data[i][1] * data[i][2];
    }

    std::cout << "One epoch of sgd on a data size " << N << " in:" << in << " out:" << out << " rows:" << rows
              << " columns:" << columns << std::endl;
    {
        boost::timer::auto_cpu_timer t;
        ex.sgd(data, label, 0.01, 32);
    }
}

/// This torture test is passed whenever it completes. It is meant to check for
/// the code stability when large number of mutations are performed
BOOST_AUTO_TEST_CASE(evaluation_speed)
{
    unsigned int N = 100;

    dcgp::kernel_set<double> kernel_set1({"sig", "tanh", "ReLu"});
    dcgp::stream(std::cout, "Function set ", kernel_set1(), "\n");
    perform_evaluations(3, 2, 2, 3, 4, 2, N, kernel_set1());
    perform_evaluations(3, 2, 10, 10, 11, 2, N, kernel_set1());
    perform_evaluations(3, 2, 20, 20, 21, 2, N, kernel_set1());
    perform_evaluations(3, 2, 1, 100, 101, 2, N, kernel_set1());
    perform_evaluations(3, 2, 2, 100, 101, 2, N, kernel_set1());
    perform_evaluations(3, 2, 3, 100, 101, 2, N, kernel_set1());

    dcgp::kernel_set<double> kernel_set2({"sig"});
    dcgp::stream(std::cout, "\nFunction set ", kernel_set2(), "\n");
    perform_evaluations(3, 2, 2, 3, 4, 2, N, kernel_set2());
    perform_evaluations(3, 2, 10, 10, 11, 2, N, kernel_set2());
    perform_evaluations(3, 2, 20, 20, 21, 2, N, kernel_set2());
    perform_evaluations(3, 2, 1, 100, 101, 2, N, kernel_set2());
    perform_evaluations(3, 2, 2, 100, 101, 2, N, kernel_set2());
    perform_evaluations(3, 2, 3, 100, 101, 2, N, kernel_set2());
}