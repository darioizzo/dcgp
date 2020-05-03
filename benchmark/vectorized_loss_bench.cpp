#include <audi/audi.hpp>
#include <boost/timer/timer.hpp>
#include <iostream>

#include <dcgp/expression.hpp>
#include <dcgp/kernel_set.hpp>

using gdual_d = audi::gdual_d;
using gdual_v = audi::gdual_v;

std::vector<std::vector<double>> transpose(const std::vector<std::vector<double>> &data)
{
    // this assumes that all inner vectors have the same size and
    // allocates space for the complete result in advance
    std::vector<std::vector<double>> result(data[0].size(), std::vector<double>(data.size()));
    for (std::vector<double>::size_type i = 0; i < data[0].size(); i++)
        for (std::vector<int>::size_type j = 0; j < data.size(); j++) {
            result[i][j] = data[j][i];
        }
    return std::move(result);
}

std::vector<std::vector<double>> random_data(unsigned N, unsigned dim, unsigned seed)
{
    // Random numbers engine
    std::default_random_engine re(seed);
    // We create the input data upfront and we do not time it.
    std::vector<double> dumb(dim);
    std::vector<std::vector<double>> retval(N, dumb);

    for (auto j = 0u; j < N; ++j) {
        for (auto i = 0u; i < dim; ++i) {
            retval[j][i] = std::uniform_real_distribution<double>(-1, 1)(re);
        }
    }
    return std::move(retval);
}

std::vector<std::vector<gdual_d>> random_data_to_gdual_d(const std::vector<std::vector<double>> &points)
{
    std::vector<gdual_d> dumb(points[0].size());
    std::vector<std::vector<gdual_d>> retval(points.size(), dumb);
    for (auto j = 0u; j < retval.size(); ++j) {
        for (auto i = 0u; i < dumb.size(); ++i) {
            retval[j][i] = gdual_d(points[j][i]);
        }
    }
    return std::move(retval);
}

std::vector<gdual_v> random_data_to_gdual_v(const std::vector<std::vector<double>> &points)
{
    std::vector<gdual_v> retval(points[0].size());
    auto pointsT = transpose(points);
    for (decltype(pointsT.size()) i = 0u; i < pointsT.size(); ++i) {
        retval[i] = gdual_v(pointsT[i]);
    }
    return std::move(retval);
}

void core_comparison(unsigned in, unsigned out, unsigned rows, unsigned columns, unsigned levels_back, unsigned arity,
                     unsigned n_eph, unsigned N, std::vector<std::string> ks)
{
    auto seed = 234234234231u;
    std::default_random_engine re(seed);
    // 1 - Instatiate a random expression in gdual_d
    unsigned order = 1u;
    dcgp::kernel_set<gdual_d> kernel_set_d(ks);
    dcgp::expression<gdual_d> ex_d(in, out, rows, columns, levels_back, arity, kernel_set_d(), n_eph, seed);
    // 2 - Set the ephemeral constants (the values by default will be 1,2,3,.... and the order zero, hence we need to
    // change)
    auto eph_val_d = ex_d.get_eph_val();
    auto eph_symb_d = ex_d.get_eph_symb();
    for (decltype(eph_val_d.size()) i = 0u; i < eph_val_d.size(); ++i) {
        eph_val_d[i] = gdual_d(i + 1u, eph_symb_d[i], order);
    }
    ex_d.set_eph_val(eph_val_d);
    // 3 - We create the input data.
    auto points_d = random_data_to_gdual_d(random_data(N, in, re()));
    auto labels_d = random_data_to_gdual_d(random_data(N, out, re()));
    // 4 - We time and compute.
    std::cout << "[gdual_d] - Performing " << N << " evaluations, in:" << in << " out:" << out << " rows:" << rows
              << " columns:" << columns << std::endl;
    {
        boost::timer::auto_cpu_timer t;
        auto res_gdual_d = ex_d.loss(points_d, labels_d, "MSE", false);
        //std::cout << "loss in gduals: " << res_gdual_d << std::endl;
        //std::cout << "expression: " << ex_d(std::vector<std::string>{"x"})[0] << std::endl;
    }

    // A - Instatiate a random expression in gdual_vd
    dcgp::kernel_set<gdual_v> kernel_set_v(ks);
    dcgp::expression<gdual_v> ex_v(in, out, rows, columns, levels_back, arity, kernel_set_v(), n_eph, seed);
    // B - Set the ephemeral constants (the values by default will be 1,2,3,.... and the order zero, hence we need to
    // change)
    auto eph_val_v = ex_v.get_eph_val();
    auto eph_symb_v = ex_v.get_eph_symb();
    for (decltype(eph_val_v.size()) i = 0u; i < eph_val_v.size(); ++i) {
        eph_val_v[i] = gdual_v(i + 1u, eph_symb_v[i], order);
    }
    ex_v.set_eph_val(eph_val_v);
    // C - We create the input data. Note that points and labels will be one only gdual_v with vector coefficients.
    std::vector<gdual_v> points_v = random_data_to_gdual_v(random_data(N, in, re()));
    std::vector<gdual_v> labels_v = random_data_to_gdual_v(random_data(N, out, re()));
    // D - We time and compute.
    std::cout << "[gdual_v] - Performing " << N << " evaluations, in:" << in << " out:" << out << " rows:" << rows
              << " columns:" << columns << std::endl;
    {
        boost::timer::auto_cpu_timer t;
        auto res_gdual_v = ex_v.loss(points_v, labels_v, dcgp::expression<gdual_v>::loss_type::MSE);
        //std::cout << "loss in gduals_v: " << res_gdual_v << std::endl;
        //std::cout << "expression: " << ex_v(std::vector<std::string>{"x"})[0] << std::endl;
        auto v = res_gdual_v.constant_cf();
        auto dvdc1 = res_gdual_v.get_derivative({1u, 0u});
        auto dvdc2 = res_gdual_v.get_derivative({0u, 1u});
        //std::cout << "loss contraction: " << std::accumulate(v.begin(), v.end(), 0.) / v.size() << std::endl;
        //std::cout << "gradient contraction (c1): " << std::accumulate(dvdc1.begin(), dvdc1.end(), 0.) / dvdc1.size()
        //          << std::endl;
        //std::cout << "gradient contraction (c2): " << std::accumulate(dvdc2.begin(), dvdc2.end(), 0.) / dvdc2.size()
        //          << std::endl;
    }
}

// This benchmark computes the loss and its gradient over a large amount of points comparing the vecotrized gdual and
// the gdual
int main()
{
    std::vector<std::string> ks{"sum", "diff", "mul", "div"};
    core_comparison(1u, 1u, 1u, 20, 21u, 2u, 2u, 10000u, ks);
    core_comparison(1u, 1u, 1u, 20, 21u, 2u, 2u, 1000u, ks);
    core_comparison(1u, 1u, 1u, 20, 21u, 2u, 2u, 100u, ks);
    core_comparison(1u, 1u, 1u, 20, 21u, 2u, 2u, 10u, ks);
    core_comparison(1u, 1u, 1u, 20, 21u, 2u, 2u, 5u, ks);

    return 0;
}
