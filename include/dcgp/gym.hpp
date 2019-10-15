#include <boost/math/constants/constants.hpp>
#include <cmath>
#include <random>
#include <vector>

namespace dcgp
{
namespace gym
{
namespace detail
{
constexpr auto pi = boost::math::constants::pi<double>();
constexpr auto e = boost::math::constants::e<double>();

typedef std::function<double(const std::vector<double> &x)> multivariate;

// Standard Problems
multivariate koza_quintic
    = [](const std::vector<double> &x) { return std::pow(x[0], 5) - 2 * std::pow(x[0], 3) + x[0]; };

// From:
// Izzo, D., Biscani, F., & Mereta, A. (2017, April).
// Differentiable genetic programming.
// In European Conference on Genetic Programming (pp. 35-51). Springer, Cham.
multivariate P1 = [](const std::vector<double> &x) { return std::pow(x[0], 5) - 2 * std::pow(x[0], 3) + x[0]; };
multivariate P2
    = [](const std::vector<double> &x) { return std::pow(x[0], 5) - pi * std::pow(x[0], 3) + 2. * pi / x[0]; };
multivariate P3
    = [](const std::vector<double> &x) { return (e * std::pow(x[0], 5) + std::pow(x[0], 3)) / (x[0] + 1.); };
multivariate P4 = [](const std::vector<double> &x) { return sin(pi * x[0]) + 1. / x[0]; };
multivariate P5 = [](const std::vector<double> &x) { return e * std::pow(x[0], 5) - pi * std::pow(x[0], 3) + x[0]; };
multivariate P6 = [](const std::vector<double> &x) { return (e * x[0] * x[0] - 1) / (pi * (x[0] + 2)); };
multivariate P7 = [](const std::vector<double> &x) { return std::cos(pi * x[0]) + std::sin(e * x[0]); };

// From:
// Vladislavleva, Ekaterina J., Guido F. Smits, and Dick Den Hertog.
// "Order of nonlinearity as a complexity measure for models generated by symbolic regression via pareto genetic
// programming." IEEE Transactions on Evolutionary Computation 13.2 (2008): 333-349. Generates data to test symbolic
// regression on 1D input-output cases.
multivariate vladi1 = [](const std::vector<double> &x) {
    return std::exp(-(x[0] - 1.) * (x[0] - 1.) / (1.2 + (x[1] - 2.5) * (x[1] - 2.5)));
};
multivariate vladi2 = [](const std::vector<double> &x) {
    return std::exp(-x[0]) * x[0] * x[0] * x[0] * std::cos(x[0]) * std::sin(x[0])
           * (std::cos(x[0]) * std::sin(x[0]) * std::sin(x[0]) - 1.);
};
multivariate vladi3 = [](const std::vector<double> &x) { return vladi2(x) * (x[1] - 5.); };
multivariate vladi4 = [](const std::vector<double> &x) {
    return 10.
           / (5.
              + std::pow((x[0] - 3.), 2) * std::pow((x[1] - 3.), 2) * std::pow((x[2] - 3.), 2) * std::pow((x[3] - 3.), 2)
                    * std::pow((x[4] - 3.), 2));
};
multivariate vladi5
    = [](const std::vector<double> &x) { return 30. * (x[0] - 1.) * (x[2] - 1.) / (x[1] * x[1] * (x[0] - 10.)); };
multivariate vladi6 = [](const std::vector<double> &x) { return 6. * std::cos(x[0] * std::sin(x[1])); };
multivariate vladi7
    = [](const std::vector<double> &x) { return (x[0] - 3.) * (x[1] - 3.) + 2 * std::sin((x[0] - 4.) * (x[1] - 4.)); };
multivariate vladi8 = [](const std::vector<double> &x) {
    return (std::pow(x[0] - 3., 4) + std::pow(x[1] - 3., 3) - (x[1] - 3.)) / (std::pow(x[1] - 2, 4.) + 10.);
};

void generate_1Ddata(std::vector<std::vector<double>> &points, std::vector<std::vector<double>> &labels, multivariate f,
                     double lb = -1, double ub = 1, unsigned N = 10)
{
    points.clear();
    labels.clear();
    for (unsigned i = 0u; i < N; ++i) {
        double x = lb + (i * (ub - lb)) / (N - 1);
        points.push_back({x});
        labels.push_back({f({x})});
    }
}
} // namespace detail
void generate_koza_quintic(std::vector<std::vector<double>> &points, std::vector<std::vector<double>> &labels)
{
    gym::detail::generate_1Ddata(points, labels, detail::koza_quintic, -1, 1, 10);
}
void generate_P1(std::vector<std::vector<double>> &points, std::vector<std::vector<double>> &labels)
{
    gym::detail::generate_1Ddata(points, labels, detail::P1, 1., 3., 10);
}
void generate_P2(std::vector<std::vector<double>> &points, std::vector<std::vector<double>> &labels)
{
    gym::detail::generate_1Ddata(points, labels, detail::P2, 0.1, 5., 10);
}
void generate_P3(std::vector<std::vector<double>> &points, std::vector<std::vector<double>> &labels)
{
    gym::detail::generate_1Ddata(points, labels, detail::P3, -0.9, 1, 10);
}
void generate_P4(std::vector<std::vector<double>> &points, std::vector<std::vector<double>> &labels)
{
    gym::detail::generate_1Ddata(points, labels, detail::P4, -1, 1, 10);
}
void generate_P5(std::vector<std::vector<double>> &points, std::vector<std::vector<double>> &labels)
{
    gym::detail::generate_1Ddata(points, labels, detail::P5, 1., 3., 10);
}
void generate_P6(std::vector<std::vector<double>> &points, std::vector<std::vector<double>> &labels)
{
    gym::detail::generate_1Ddata(points, labels, detail::P6, -2.1, 1., 10);
}
void generate_P7(std::vector<std::vector<double>> &points, std::vector<std::vector<double>> &labels)
{
    gym::detail::generate_1Ddata(points, labels, detail::P7, -1, 1, 10);
}
void generate_vladi1(std::vector<std::vector<double>> &points, std::vector<std::vector<double>> &labels)
{
    points.clear();
    labels.clear();
    std::mt19937 mt(32);
    std::uniform_real_distribution<double> dist(0.3, 4);
    for (unsigned i = 0; i < 100u; ++i) {
        std::vector<double> point = {dist(mt), dist(mt)};
        points.push_back(point);
        labels.push_back({detail::vladi1(point)});
    }
}
void generate_vladi2(std::vector<std::vector<double>> &points, std::vector<std::vector<double>> &labels)
{
    gym::detail::generate_1Ddata(points, labels, detail::vladi2, 0.05, 10., 100);
}
void generate_vladi3(std::vector<std::vector<double>> &points, std::vector<std::vector<double>> &labels)
{
    points.clear();
    labels.clear();
    std::mt19937 mt(32);
    std::uniform_real_distribution<double> dist(0.05, 10);
    for (unsigned i = 0; i < 601u; ++i) {
        std::vector<double> point = {dist(mt), dist(mt)};
        points.push_back(point);
        labels.push_back({detail::vladi3(point)});
    }
}
void generate_vladi4(std::vector<std::vector<double>> &points, std::vector<std::vector<double>> &labels)
{
    points.clear();
    labels.clear();
    std::mt19937 mt(32);
    std::uniform_real_distribution<double> dist(0.05, 6.05);
    for (unsigned i = 0; i < 1024u; ++i) {
        std::vector<double> point = {dist(mt), dist(mt), dist(mt), dist(mt), dist(mt)};
        points.push_back(point);
        labels.push_back({detail::vladi4(point)});
    }
}
void generate_vladi5(std::vector<std::vector<double>> &points, std::vector<std::vector<double>> &labels)
{
    points.clear();
    labels.clear();
    std::mt19937 mt(32);
    std::uniform_real_distribution<double> dist(0.05, 2);
    std::uniform_real_distribution<double> dist1(1, 2);
    for (unsigned i = 0; i < 300u; ++i) {
        std::vector<double> point = {dist(mt), dist1(mt), dist(mt)};
        points.push_back(point);
        labels.push_back({detail::vladi5(point)});
    }
}
void generate_vladi6(std::vector<std::vector<double>> &points, std::vector<std::vector<double>> &labels)
{
    points.clear();
    labels.clear();
    std::mt19937 mt(32);
    std::uniform_real_distribution<double> dist(0.05, 2);
    for (unsigned i = 0; i < 30u; ++i) {
        std::vector<double> point = {dist(mt), dist(mt)};
        points.push_back(point);
        labels.push_back({detail::vladi6(point)});
    }
}
void generate_vladi7(std::vector<std::vector<double>> &points, std::vector<std::vector<double>> &labels)
{
    points.clear();
    labels.clear();
    std::mt19937 mt(32);
    std::uniform_real_distribution<double> dist(0.05, 6.05);
    for (unsigned i = 0; i < 300u; ++i) {
        std::vector<double> point = {dist(mt), dist(mt)};
        points.push_back(point);
        labels.push_back({detail::vladi7(point)});
    }
}
void generate_vladi8(std::vector<std::vector<double>> &points, std::vector<std::vector<double>> &labels)
{
    points.clear();
    labels.clear();
    std::mt19937 mt(32);
    std::uniform_real_distribution<double> dist(0.05, 6.05);
    for (unsigned i = 0; i < 50u; ++i) {
        std::vector<double> point = {dist(mt), dist(mt)};
        points.push_back(point);
        labels.push_back({detail::vladi8(point)});
    }
}
} // namespace gym
} // namespace dcgp