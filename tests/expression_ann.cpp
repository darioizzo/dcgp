#include <random>
#define BOOST_TEST_MODULE dcgp_expression_ann_test
#include <algorithm>
#include <audi/back_compatibility.hpp>
#include <audi/io.hpp>
#include <boost/test/unit_test.hpp>
#include <stdexcept>

#include <dcgp/dcgp.hpp>

using namespace dcgp;

BOOST_AUTO_TEST_CASE(construction)
{
    // Random seed
    std::random_device rd;
    // Kernel functions
    kernel_set<double> ann_set({"tanh"});
    expression_ann<double> ex(1, 1, 1, 2, 1, 1, ann_set(), rd());
    // We test that all weights are set to 1 and biases to 0
    auto ws = ex.get_weights();
    auto bs = ex.get_biases();
    BOOST_CHECK(std::all_of(ws.begin(), ws.end(), [](unsigned el) { return el == 1u; }));
    BOOST_CHECK(std::all_of(bs.begin(), bs.end(), [](unsigned el) { return el == 0u; }));

    kernel_set<double> ann_set_malformed1({"tanh", "sin"});
    kernel_set<double> ann_set_malformed2({"cos", "sig"});
    kernel_set<double> ann_set_malformed3({"ReLu", "sum"});

    BOOST_CHECK_THROW((expression_ann<double>{1, 1, 1, 2, 1, 1, ann_set_malformed1(), rd()}), std::invalid_argument);
    BOOST_CHECK_THROW((expression_ann<double>{1, 1, 1, 2, 1, 1, ann_set_malformed2(), rd()}), std::invalid_argument);
    BOOST_CHECK_THROW((expression_ann<double>{1, 1, 1, 2, 1, 1, ann_set_malformed3(), rd()}), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(parenthesis)
{
    {
        // We test a simple arity 1 row 1 dCGP-ANN
        // Random seed
        std::random_device rd;
        // Kernel functions
        kernel_set<double> ann_set({"tanh"});
        expression_ann<double> ex(1, 1, 1, 2, 1, 1, ann_set(), rd());
        ex.set_weights({0.1, 0.2});
        ex.set_biases({0.3, 0.4});
        auto res = ex({0.23})[0];
        auto ground_truth = std::tanh(0.4 + 0.2 * std::tanh(0.23 * 0.1 + 0.3));
        BOOST_CHECK_CLOSE(res, ground_truth, 1e-13);
    }
    {
        // We test a simple arity 2 row 1 dCGP-ANN
        // Random seed
        std::random_device rd;
        // Kernel functions
        kernel_set<double> ann_set({"tanh"});
        expression_ann<double> ex(1, 1, 1, 2, 1, 2, ann_set(), rd());
        ex.set_weights({0.1, 0.2, 0.3, 0.4});
        ex.set_biases({0.5, 0.6});
        auto res = ex({0.23})[0];
        auto n1 = std::tanh(0.23 * 0.1 + 0.23 * 0.2 + 0.5);
        auto ground_truth = std::tanh(0.3 * n1 + 0.4 * n1 + 0.6);
        BOOST_CHECK_CLOSE(res, ground_truth, 1e-13);
    }
    {
        // We test a arity 2 row 2 column 2 dCGP-ANN
        // Random seed
        std::random_device rd;
        // Kernel functions
        kernel_set<double> ann_set({"tanh"});
        expression_ann<double> ex(1, 1, 2, 2, 1, 2, ann_set(), rd());
        ex.set_weights({0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8});
        ex.set_biases({0.9, 1.1, 1.2, 1.3});
        ex.set({0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 1, 2, 3});
        auto res = ex({0.23})[0];
        auto n0 = 0.23;
        auto n1 = std::tanh(0.1 * n0 + 0.2 * n0 + 0.9);
        auto n2 = std::tanh(0.3 * n0 + 0.4 * n0 + 1.1);
        auto ground_truth = std::tanh(0.5 * n1 + 0.6 * n2 + 1.2);
        BOOST_CHECK_CLOSE(res, ground_truth, 1e-13);
    }
}

BOOST_AUTO_TEST_CASE(mse)
{
    for (auto i = 0u; i < 1000; ++i) {

        {
            print(i, "\n");
            // We test a arity 10 row 100 column 3, 3 input 2 outputs dCGP-ANN
            // Random seed
            std::random_device rd;
            // Kernel functions
            kernel_set<double> ann_set({"sig", "tanh", "ReLu"});
            expression_ann<double> ex(3, 2, 100, 3, 1, 10, ann_set(), rd());
            ex.randomise_weights();
            ex.randomise_biases();
            auto orig_w = ex.get_weights();
            auto orig_b = ex.get_biases();
            // Numerical derivative eps (low precision but more reliable)
            auto eps = 1e-4;
            // Input value
            auto in = std::vector<double>(3, 0.22);
            // Output value desired (supervised signal)
            auto out = std::vector<double>(2, 0.23);
            // Compute mse and the gradients
            auto bp = ex.mse(in, out);
            // We check against numerical diff within 20% of accuracy
            // first the weights
            for (decltype(ex.get_weights().size()) i = 0u; i < ex.get_weights().size(); ++i) {
                ex.set_weights(orig_w);
                auto tmp = ex.get_weight(i);
                auto h = std::max(1., std::abs(tmp)) * eps;
                ex.set_weight(i, tmp + h);
                auto val = (ex(in)[0] - out[0]) * (ex(in)[0] - out[0]) + (ex(in)[1] - out[1]) * (ex(in)[1] - out[1]);
                ex.set_weight(i, tmp - h);
                auto val2 = (ex(in)[0] - out[0]) * (ex(in)[0] - out[0]) + (ex(in)[1] - out[1]) * (ex(in)[1] - out[1]);

                auto flag = std::abs(((val - val2) / 2. / h - std::get<1>(bp)[i]) / ((val - val2) / 2. / h));
                // if (std::get<2>(bp)[i] != 0) {
                //    print(flag, "\n");
                // }
                if (flag > 0.2) {
                    print("Epic Fail:\n ");
                    print(ex);
                    auto steps = std::vector<double>{
                        1.,        0.1,        0.01,        0.001,        0.0001,        0.00001,       0.000001,
                        0.0000001, 0.00000001, 0.000000001, 0.0000000001, 1.00000000001, 1.000000000001};
                    for (auto hh : steps) {
                        ex.set_weights(orig_w);
                        ex.set_weight(i, tmp + hh);
                        auto val
                            = (ex(in)[0] - out[0]) * (ex(in)[0] - out[0]) + (ex(in)[1] - out[1]) * (ex(in)[1] - out[1]);
                        ex.set_weight(i, tmp - hh);
                        auto val2
                            = (ex(in)[0] - out[0]) * (ex(in)[0] - out[0]) + (ex(in)[1] - out[1]) * (ex(in)[1] - out[1]);
                        print("Num Grad: ", (val - val2) / 2. / hh, ", h used: ", h, ", bias: ", tmp, "\n");
                    }
                }

                BOOST_CHECK_CLOSE((val - val2) / 2. / h, std::get<1>(bp)[i], 5.);
            }
            // then the biases
            ex.set_weights(orig_w);
            for (decltype(ex.get_biases().size()) i = 0u; i < ex.get_biases().size(); ++i) {
                ex.set_biases(orig_b);
                auto tmp = ex.get_bias(i);
                auto h = std::max(1., std::abs(tmp)) * eps;
                ex.set_bias(i, tmp + h);
                auto val = (ex(in)[0] - out[0]) * (ex(in)[0] - out[0]) + (ex(in)[1] - out[1]) * (ex(in)[1] - out[1]);
                ex.set_bias(i, tmp - h);
                auto val2 = (ex(in)[0] - out[0]) * (ex(in)[0] - out[0]) + (ex(in)[1] - out[1]) * (ex(in)[1] - out[1]);
                auto flag = std::abs(((val - val2) / 2. / h - std::get<2>(bp)[i]) / ((val - val2) / 2. / h));
                // if (std::get<2>(bp)[i] != 0) {
                //    print(flag, "\n");
                // }
                if (flag > 0.2) {
                    print("Epic Fail:\n ");
                    print(ex);
                    auto steps = std::vector<double>{
                        1.,        0.1,        0.01,        0.001,        0.0001,        0.00001,       0.000001,
                        0.0000001, 0.00000001, 0.000000001, 0.0000000001, 1.00000000001, 1.000000000001};
                    for (auto hh : steps) {
                        ex.set_biases(orig_b);
                        ex.set_bias(i, tmp + hh);
                        auto val
                            = (ex(in)[0] - out[0]) * (ex(in)[0] - out[0]) + (ex(in)[1] - out[1]) * (ex(in)[1] - out[1]);
                        ex.set_bias(i, tmp - hh);
                        auto val2
                            = (ex(in)[0] - out[0]) * (ex(in)[0] - out[0]) + (ex(in)[1] - out[1]) * (ex(in)[1] - out[1]);
                        print("Num Grad: ", (val - val2) / 2. / hh, ", h used: ", h, ", bias: ", tmp, "\n");
                    }
                }

                BOOST_CHECK_CLOSE((val - val2) / 2 / h, std::get<2>(bp)[i], 5.);
            }
        }
    }
}
