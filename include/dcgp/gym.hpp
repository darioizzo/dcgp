#include <boost/math/constants/constants.hpp>
#include <cmath>
#include <vector>

namespace dcgp
{
namespace gym
{
namespace detail
{
constexpr auto pi = boost::math::constants::pi<double>();
constexpr auto e = boost::math::constants::e<double>();

// Standard Problems
std::function<double(double)> koza_quintic_fun = [](double x) { return std::pow(x, 5) - 2 * std::pow(x, 3) + x; };
// From our paper on dCGP
std::function<double(double)> P1_fun = [](double x) { return std::pow(x, 5) - 2 * std::pow(x, 3) + x; };
std::function<double(double)> P2_fun = [](double x) { return std::pow(x, 5) - pi * std::pow(x, 3) + 2 * pi / x; };
std::function<double(double)> P3_fun = [](double x) { return (e * std::pow(x, 5) + std::pow(x, 3)) / (x + 1); };
std::function<double(double)> P4_fun = [](double x) { return sin(pi * x) + 1 / x; };
std::function<double(double)> P5_fun = [](double x) { return e * std::pow(x, 5) - pi * std::pow(x, 3) + x; };
std::function<double(double)> P6_fun = [](double x) { return (e * x * x - 1) / (pi * (x + 2)); };
std::function<double(double)> P7_fun = [](double x) { return std::cos(pi * x) + std::sin(e * x); };

// Generates data to test symbolic regression on 1D input-output cases
void generate_1Ddata(std::vector<std::vector<double>> &points, std::vector<std::vector<double>> &labels,
                     std::function<double(double)> f, double lb = -1, double ub = 1, unsigned N = 10)
{
    points.clear();
    labels.clear();
    for (unsigned i = 0u; i < N; ++i) {
        double x = lb + i / (N-1) * (ub - lb);
        points.push_back({x});
        labels.push_back({f(x)});
    }
}
} // namespace detail
void generate_koza_quintic(std::vector<std::vector<double>> &points, std::vector<std::vector<double>> &labels)
{
    gym::detail::generate_1Ddata(points, labels, detail::koza_quintic_fun, -1, 1, 10);
}
void generate_P1(std::vector<std::vector<double>> &points, std::vector<std::vector<double>> &labels)
{
    gym::detail::generate_1Ddata(points, labels, detail::P1_fun, 1., 3., 10);
}
void generate_P2(std::vector<std::vector<double>> &points, std::vector<std::vector<double>> &labels)
{
    gym::detail::generate_1Ddata(points, labels, detail::P2_fun, 0.1, 5., 10);
}
void generate_P3(std::vector<std::vector<double>> &points, std::vector<std::vector<double>> &labels)
{
    gym::detail::generate_1Ddata(points, labels, detail::P3_fun, -0.9, 1, 10);
}
void generate_P4(std::vector<std::vector<double>> &points, std::vector<std::vector<double>> &labels)
{
    gym::detail::generate_1Ddata(points, labels, detail::P4_fun, -1, 1, 10);
}
void generate_P5(std::vector<std::vector<double>> &points, std::vector<std::vector<double>> &labels)
{
    gym::detail::generate_1Ddata(points, labels, detail::P5_fun, 1., 3., 10);
}
void generate_P6(std::vector<std::vector<double>> &points, std::vector<std::vector<double>> &labels)
{
    gym::detail::generate_1Ddata(points, labels, detail::P6_fun, -2.1, 1., 10);
}
void generate_P7(std::vector<std::vector<double>> &points, std::vector<std::vector<double>> &labels)
{
    gym::detail::generate_1Ddata(points, labels, detail::P7_fun, -1, 1, 10);
}
} // namespace gym
} // namespace dcgp