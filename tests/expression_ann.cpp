#include <random>
#define BOOST_TEST_MODULE dcgp_expression_ann_test
#include <algorithm>
#include <audi/back_compatibility.hpp>
#include <audi/io.hpp>
#include <boost/test/unit_test.hpp>
#include <stdexcept>
#include <boost/lexical_cast.hpp>

#include <dcgp/dcgp.hpp>

using namespace dcgp;
/*
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

BOOST_AUTO_TEST_CASE(sgd)
{
    // Random numbers stuff
    std::random_device rd;
    std::mt19937 gen{rd()};
    std::normal_distribution<> norm(0.,1.);

    // Kernel functions
    kernel_set<double> ann_set({"sig", "tanh", "ReLu"});
    expression_ann<double> ex(3, 2, 100, 3, 1, 10, ann_set(), rd());
    ex.randomise_weights();
    ex.randomise_biases();
    std::vector<std::vector<double>> data(100, {0.,0.,0.});
    std::vector<std::vector<double>> label(100, {0., 0.});
    for (auto &item : data) {
        std::generate(item.begin(), item.end(), [&norm, &gen](){return norm(gen);});
    }
    for (auto i = 0u; i < label.size(); ++i) {
        label[i][0] = 1./5.*std::cos(data[i][0]+data[i][1]+data[i][2]) - data[i][0]*data[i][1];
        label[i][1] = data[i][0]*data[i][1]*data[i][2];
    }
    double tmp = 0.;
    for (auto i = 0u; i < data.size(); ++i) {
        tmp += std::get<0>(ex.mse(data[i], label[i]));
    }
    tmp /= static_cast<double>(data.size());
    print("Start: ", tmp, "\n");
    print("Start: ", std::get<0>(ex.mse(data,label)), "\n");
    for (auto j = 0u; j < 10; ++j) {
        ex.sgd(data, label, 0.1, 32);
        tmp = 0.;
        for (auto i = 0u; i < data.size(); ++i) {
            tmp += std::get<0>(ex.mse(data[i], label[i]));
        }
        tmp /= static_cast<double>(data.size());
        print("Then (", j, "): ", tmp, "\n");
    }
}
*/

BOOST_AUTO_TEST_CASE(mse)
{
    {
        // Random seed
        std::random_device rd;
        // Random generator
        auto seed = rd();
        //seed = 1186804711u;
        print("seed: ", seed, '\n');
        std::mt19937 gen(seed);
        // Random distributions
        std::uniform_int_distribution<> in_size(1, 5);
        std::uniform_int_distribution<> out_size(1, 5);
        std::uniform_int_distribution<> row(1, 20);
        std::uniform_int_distribution<> col(1, 20);
        std::uniform_int_distribution<> lb(1, 20);
        std::uniform_int_distribution<> ar(2, 10);
        std::normal_distribution<> norm{0., 1.};
        std::uniform_int_distribution<unsigned> random_seed(2, 1654636360u);

        // Kernel functions
        kernel_set<double> ann_set({"sig", "tanh", "ReLu"});
        //expression_ann<double> ex(3, 2, 100, 3, 2, 10, ann_set(), rd());
        // the dCGPANN
        expression_ann<double> ex(in_size(gen), out_size(gen), row(gen), col(gen), lb(gen), ar(gen), ann_set(),random_seed(gen));
        //expression_ann<double> ex(2, 2, 2, 2, 2, 2, ann_set(), random_seed(gen));
        // Since weights and biases are, by default, set to ones we randomize them
        ex.randomise_weights(0, 1., random_seed(gen));
        ex.randomise_biases(0, 1., random_seed(gen));

        auto orig_w = ex.get_weights();
        auto orig_b = ex.get_biases();
        // Numerical derivative eps (low precision but more reliable)
        auto eps = 1e-4;
        // Input value
        auto in = std::vector<double>(ex.get_n(), norm(gen));
        // Output value desired (supervised signal)
        auto out = std::vector<double>(ex.get_m(), norm(gen));
        // Compute mse and the gradients
        auto bp = ex.mse(in, out);
        // We check against numerical diff within accuracy
        print("In size: ", ex.get_n(), "\nOut size: ", ex.get_m(), "\nRows: ", ex.get_rows(), "\nCols: ", ex.get_cols(),
              "\nLevel-backs: ", ex.get_levels_back(), "\nArity: ", ex.get_arity(), "\n");
        print("\nNumber of Weights: ", orig_w.size(), ", number of biases: ", orig_b.size(), "\n");
        print("\nChecking Weights\n");
        // first the weights
        ex.set_weights(orig_w);
        ex.set_biases(orig_b);
        for (decltype(ex.get_weights().size()) i = 0u; i < ex.get_weights().size(); ++i) {

        //for (decltype(ex.get_weights().size()) i = 2u; i < 3; ++i) {

        
            ex.set_weights(orig_w);
            auto tmp = ex.get_weight(i);
            auto h = std::max(1., std::abs(tmp)) * eps;
            ex.set_weight(i, tmp + h);
            auto res = ex(in);
            std::transform(res.begin(), res.end(), out.begin(), res.begin(), [](double a, double b) { return a - b; });
            auto val = std::inner_product(res.begin(), res.end(), res.begin(), 0.0);
            ex.set_weight(i, tmp - h);
            res = ex(in);
            std::transform(res.begin(), res.end(), out.begin(), res.begin(), [](double a, double b) { return a - b; });
            auto val2 = std::inner_product(res.begin(), res.end(), res.begin(), 0.0);

            // Since numerical differentiation sucks we look (brute force) for a better step for numerical
            // differentiation and hope for the best if the default h failed
            auto abs_diff = std::abs(((val - val2) / 2. / h - std::get<1>(bp)[i]));
            auto rel_diff = abs_diff / std::abs(std::get<1>(bp)[i]);
            auto best = rel_diff;
            auto bval = val;
            auto bval2 = val2;
            std::string debug_stream;
            if (rel_diff > 0.05) {
                h = 100.;
                for (auto j = 0u; j < 12; ++j) {
                    ex.set_weights(orig_w);
                    tmp = ex.get_weight(i);
                    h = h * 0.1; // will generate 10, 1, 0.1, 0.001, ...., 0.000000001
                    ex.set_weight(i, tmp + h);
                    res = ex(in);
                    std::transform(res.begin(), res.end(), out.begin(), res.begin(),
                                   [](double a, double b) { return a - b; });
                    val = std::inner_product(res.begin(), res.end(), res.begin(), 0.0);
                    ex.set_weight(i, tmp - h);
                    res = ex(in);
                    std::transform(res.begin(), res.end(), out.begin(), res.begin(),
                                   [](double a, double b) { return a - b; });
                    val2 = std::inner_product(res.begin(), res.end(), res.begin(), 0.0);
                    abs_diff = std::abs(((val - val2) / 2. / h - std::get<1>(bp)[i]));
                    rel_diff = abs_diff / std::abs(std::get<1>(bp)[i]);
                    debug_stream += "\nh: " + boost::lexical_cast<std::string>(h) + " rel_diff: " + boost::lexical_cast<std::string>(rel_diff) + "\n";
                    debug_stream += "h: " + boost::lexical_cast<std::string>(h) + " abs_diff: " + boost::lexical_cast<std::string>(abs_diff) + "\n";
                    debug_stream += "val: " + boost::lexical_cast<std::string>(val) + " val2: " + boost::lexical_cast<std::string>(val2) + "\n";
                    debug_stream += "der: " + boost::lexical_cast<std::string>((val - val2) / 2. / h) + "\n";
                    debug_stream += "ana der: " +  boost::lexical_cast<std::string>(std::get<1>(bp)[i]) + "\n";
                    if (rel_diff < best) {
                        best = rel_diff;
                        bval = val;
                        bval2 = val2;
                    }
                    if (rel_diff < 0.05 || abs_diff == 0) {
                        break;
                    }
                }
                debug_stream += "Weight id: " + boost::lexical_cast<std::string>(i) + "\n";
                debug_stream += "Difference fixed to: " + boost::lexical_cast<std::string>(best) + "\n";
            }
            if (bval != bval2) {
                if (best < 0.05 || abs_diff < 1e-8) {
                    BOOST_CHECK(best < 0.05 || abs_diff < 1e-8);
                } else {
                    print(debug_stream, "\n");
                    print("\nExpression: ", ex(std::vector<std::string>{"x", "y"}));
                    print("\nExpression: ", ex(in));
                    print("\nChromosome: [");
                    for (auto i : ex.get()) {
                        print(i, ",");
                    }
                    print("]\n");
                    print("Active nodes: [");
                    for (auto i : ex.get_active_nodes()) {
                        print(i, ",");
                    }
                    print("]\n");
                    print("\nWeights: [");
                    for (auto i : ex.get_weights()) {
                        print(i, ",");
                    }
                    print("]");
                    print("\nBiases: [");
                    for (auto i : ex.get_biases()) {
                        print(i, ",");
                    }
                    print("]");
                    print("\nWeights Derivatives: [");
                    for (auto i : std::get<1>(bp)) {
                        print(i, ",");
                    }
                    print("]\n");
                    print("Biases Derivatives: [");
                    for (auto i : std::get<2>(bp)) {
                        print(i, ",");
                    }
                    print("]\n\n");
                    print("Input: [");
                    for (auto i : in) {
                        print(i, ",");
                    }
                    print("]\n");
                    print("Output desired: [");
                    for (auto i : out) {
                        print(i, ",");
                    }
                    print("]\n");
                    print("bval ", bval, "\n");
                    print("bval2 ", bval2, "\n");
                    print("Equal?: ", bval==bval2, "\n");
                    throw std::invalid_argument("AAARGHHH");
                }
            } else {
                BOOST_CHECK(std::abs(std::get<1>(bp)[i]) < 1e-8);
            }
        }
        print("Checking Biases\n");
        // then the biases
        ex.set_weights(orig_w);
        ex.set_biases(orig_b);
        for (decltype(ex.get_biases().size()) i = 0u; i < ex.get_biases().size(); ++i) {
            ex.set_biases(orig_b);
            auto tmp = ex.get_bias(i);
            auto h = std::max(1., std::abs(tmp)) * eps;
            ex.set_bias(i, tmp + h);
            auto res = ex(in);
            std::transform(res.begin(), res.end(), out.begin(), res.begin(), [](double a, double b) { return a - b; });
            auto val = std::inner_product(res.begin(), res.end(), res.begin(), 0.0);
            ex.set_bias(i, tmp - h);
            res = ex(in);
            std::transform(res.begin(), res.end(), out.begin(), res.begin(), [](double a, double b) { return a - b; });
            auto val2 = std::inner_product(res.begin(), res.end(), res.begin(), 0.0);
            // BOOST_CHECK_CLOSE((val - val2) / 2 / h, std::get<2>(bp)[i], 5.);
        }
    }
}
