#include <random>
#define BOOST_TEST_MODULE symbolic_regresion_perf_test
#include <algorithm>
#include <audi/back_compatibility.hpp>
#include <audi/io.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/timer/timer.hpp>

#include <dcgp/gym.hpp>
#include <dcgp/kernel_set.hpp>
#include <dcgp/problems/symbolic_regression.hpp>

using namespace dcgp;

BOOST_AUTO_TEST_CASE(fitness_gradient_speed)
{
    kernel_set<double> basic_set({"sum", "diff", "mul", "div"});
    std::vector<std::vector<double>> points, labels;
    gym::generate_P4(points, labels);
    symbolic_regression udp(points, labels, 2, 100, 30, 2, basic_set(), 5u, 0u);
    pagmo::population pop(udp, 100, 32u);

    // We measure the speed on N fitness evaluations
    pagmo::print("Fitness evaluation speed on 100 calls: \n");
    {
        boost::timer::auto_cpu_timer t;
        for (decltype(pop.size()) i = 0u; i < pop.size(); ++i) {
            udp.fitness(pop.get_x()[i]);
        }
    }
    // We measure the speed on N gradient evaluations
    pagmo::print("Gradient evaluation speed on 100 calls: \n");
    {
        boost::timer::auto_cpu_timer t;
        for (decltype(pop.size()) i = 0u; i < pop.size(); ++i) {
            udp.gradient(pop.get_x()[i]);
        }
    }
}

BOOST_AUTO_TEST_CASE(pretty_speed)
{
    kernel_set<double> basic_set({"sum", "diff", "mul", "div"});
    std::vector<std::vector<double>> points, labels;
    gym::generate_P4(points, labels);
    symbolic_regression udp(points, labels, 2, 100, 30, 2, basic_set(), 5u, 0u);
    pagmo::population pop(udp, 100, 32u);

    // We measure the speed of pretty
    pagmo::print("\n\nPretty called on 100 chromosomes: \n");
    {
        boost::timer::auto_cpu_timer t;
        for (decltype(pop.size()) i = 0u; i < pop.size(); ++i) {
            udp.pretty(pop.get_x()[i]);
        }
    }
    // We measure the speed on N gradient evaluations
    pagmo::print("Prettier called on 100 chromosomes: \n");
    {
        boost::timer::auto_cpu_timer t;
        for (decltype(pop.size()) i = 0u; i < pop.size(); ++i) {
            udp.prettier(pop.get_x()[i]);
        }
    }
    double anomalies = 0;
    double diff = 0;
    // We measure the length difference
    pagmo::print("\nWe measure the difference in length of the expressions: \n");
    for (decltype(pop.size()) i = 0u; i < pop.size(); ++i) {
        auto s = udp.prettier(pop.get_x()[i]);
        auto l = udp.pretty(pop.get_x()[i]);
        diff += (static_cast<double>(l.size()) - static_cast<double>(s.size()));
        if (l.size() < s.size()) {
            anomalies+=1;
        }
    }
    pagmo::print("Average length difference: ", diff / 100., "\n", "Number of anomalies (prettier less pretty than pretty): ", anomalies, "\n");
}
