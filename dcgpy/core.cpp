#include <audi/audi.hpp>
#include <boost/python.hpp>
#include <boost/python/tuple.hpp>
#include <functional> //std::function
#include <sstream>
#include <string>
#include <vector>

#include <dcgp/expression.hpp>
#include <dcgp/expression_ann.hpp>
#include <dcgp/expression_weighted.hpp>
#include <dcgp/gym.hpp>
#include <dcgp/kernel.hpp>
#include <dcgp/kernel_set.hpp>

// See: https://docs.scipy.org/doc/numpy/reference/c-api.array.html#importing-the-api
// In every cpp file We need to make sure this is included before everything else,
// with the correct #defines.
#define PY_ARRAY_UNIQUE_SYMBOL dcgpy_ARRAY_API

#include "common_utils.hpp"
#include "docstrings.hpp"

using namespace dcgp;
using namespace dcgpy;
using namespace audi;
namespace bp = boost::python;

template <typename T>
void expose_kernel(const std::string &type)
{
    std::string class_name = "kernel_" + type;
    bp::class_<kernel<T>>(class_name.c_str(), "The function defining the generic CGP node", bp::no_init)
        .def("__init__",
             bp::make_constructor(
                 +[](const bp::object &obj1, const bp::object &obj2, const std::string &name) {
                     std::function<T(const std::vector<T> &)> my_function = [obj1](const std::vector<T> &x) {
                         T in = bp::extract<T>(obj1(v_to_l(x)));
                         return in;
                     };
                     std::function<std::string(const std::vector<std::string> &)> my_print_function
                         = [obj2](const std::vector<std::string> &x) {
                               std::string in = bp::extract<std::string>(obj2(v_to_l(x)));
                               return in;
                           };
                     return ::new kernel<T>(my_function, my_print_function, name);
                 },
                 bp::default_call_policies(), (bp::arg("callable_f"), bp::arg("callable_s"), bp::arg("name"))),
             kernel_init_doc(type).c_str())
        .def(
            "__call__",
            +[](kernel<T> &instance, const bp::object &in) {
                try {
                    auto v = l_to_v<T>(in);
                    return bp::object(instance(v));
                } catch (...) {
                    PyErr_Clear();
                    auto v = l_to_v<std::string>(in);
                    return bp::object(instance(v));
                }
            })
        .def(
            "__repr__", +[](const kernel<T> &instance) -> std::string {
                std::ostringstream oss;
                oss << instance;
                return oss.str();
            });
    ;
}

template <typename T>
kernel<T> wrap_operator(const kernel_set<T> &ks, typename std::vector<dcgp::kernel<T>>::size_type idx)
{
    return ks[idx];
}

template <typename T>
void expose_kernel_set(std::string type)
{
    std::string class_name = "kernel_set_" + type;
    bp::class_<kernel_set<T>>(class_name.c_str(),
                              "Helper to construct a set of kernel functions from their common name", bp::no_init)
        .def("__init__",
             bp::make_constructor(
                 +[](const bp::object &obj1) {
                     auto a = l_to_v<std::string>(obj1);
                     return ::new kernel_set<T>(a);
                 },
                 bp::default_call_policies(), (bp::arg("kernels"))),
             kernel_set_init_doc(type).c_str())
        .def(
            "__call__", +[](kernel_set<T> &instance) { return v_to_l(instance()); })
        .def(
            "__repr__",
            +[](const kernel_set<T> &instance) -> std::string {
                std::ostringstream oss;
                oss << instance;
                return oss.str();
            })
        .def("push_back", (void (kernel_set<T>::*)(std::string)) & kernel_set<T>::push_back,
             kernel_set_push_back_str_doc().c_str(), bp::arg("kernel_name"))
        .def("push_back", (void (kernel_set<T>::*)(const kernel<T> &)) & kernel_set<T>::push_back,
             kernel_set_push_back_ker_doc(type).c_str(), bp::arg("kernel"))
        .def("__getitem__", &wrap_operator<T>);
}

template <typename T>
void expose_expression(std::string type)
{
    std::string class_name = "expression_" + type;
    bp::class_<expression<T>>(class_name.c_str(), "A CGP expression", bp::no_init)
        // Constructor with seed
        .def("__init__",
             bp::make_constructor(
                 +[](unsigned in, unsigned out, unsigned rows, unsigned cols, unsigned levelsback,
                     const bp::object &arity, const bp::object &kernels, unsigned n_eph, unsigned seed) {
                     auto kernels_v = l_to_v<kernel<T>>(kernels);
                     bp::extract<unsigned> is_int(arity);
                     if (is_int.check()) { // arity is passed as an integer
                         unsigned ar = bp::extract<unsigned>(arity);
                         return ::new expression<T>(in, out, rows, cols, levelsback, ar, kernels_v, n_eph, seed);
                     } else { // arity is passed as something else, a list is assumed
                         auto varity = l_to_v<unsigned>(arity);
                         return ::new expression<T>(in, out, rows, cols, levelsback, varity, kernels_v, n_eph, seed);
                     }
                 },
                 bp::default_call_policies(),
                 (bp::arg("inputs"), bp::arg("outputs"), bp::arg("rows"), bp::arg("cols"), bp::arg("levels_back"),
                  bp::arg("arity"), bp::arg("kernels"), bp::arg("n_eph"), bp::arg("seed"))),
             expression_init_doc(type).c_str())
        // Constructor with no seed
        .def("__init__",
             bp::make_constructor(
                 +[](unsigned in, unsigned out, unsigned rows, unsigned cols, unsigned levelsback,
                     const bp::object &arity, const bp::object &kernels, unsigned n_eph) {
                     auto kernels_v = l_to_v<kernel<T>>(kernels);
                     bp::extract<unsigned> is_int(arity);
                     if (is_int.check()) { // arity is passed as an integer
                         unsigned ar = bp::extract<unsigned>(arity);
                         return ::new expression<T>(in, out, rows, cols, levelsback, ar, kernels_v, n_eph,
                                                    std::random_device()());
                     } else { // arity is passed as something else, a list is assumed
                         auto varity = l_to_v<unsigned>(arity);
                         return ::new expression<T>(in, out, rows, cols, levelsback, varity, kernels_v, n_eph,
                                                    std::random_device()());
                     }
                 },
                 bp::default_call_policies(),
                 (bp::arg("inputs"), bp::arg("outputs"), bp::arg("rows"), bp::arg("cols"), bp::arg("levels_back"),
                  bp::arg("arity"), bp::arg("kernels"), bp::arg("n_eph"))),
             expression_init_doc(type).c_str())
        .def(
            "__repr__",
            +[](const expression<T> &instance) -> std::string {
                std::ostringstream oss;
                oss << instance;
                return oss.str();
            })
        .def(
            "__call__",
            +[](const expression<T> &instance, const bp::object &in) {
                try {
                    auto v = l_to_v<T>(in);
                    return v_to_l(instance(v));
                } catch (...) {
                    PyErr_Clear();
                    auto v = l_to_v<std::string>(in);
                    return v_to_l(instance(v));
                }
            })
        .def(
            "set", +[](expression<T> &instance, const bp::object &in) { instance.set(l_to_v<unsigned>(in)); },
            expression_set_doc().c_str(), bp::arg("chromosome"))
        .def("set_f_gene", &expression<T>::set_f_gene, expression_set_f_gene_doc().c_str(),
             (bp::arg("node_id"), bp::arg("f_id")))
        .def(
            "get", +[](const expression<T> &instance) { return v_to_l(instance.get()); },
            "Gets the expression chromosome")
        .def(
            "get_lb", +[](const expression<T> &instance) { return v_to_l(instance.get_lb()); },
            "Gets the lower bounds of the chromosome")
        .def(
            "get_ub", +[](const expression<T> &instance) { return v_to_l(instance.get_ub()); },
            "Gets the upper bounds of the chromosome")
        .def(
            "get_active_genes", +[](const expression<T> &instance) { return v_to_l(instance.get_active_genes()); },
            "Gets the idx of the active genes in the current chromosome (numbering is from 0)")
        .def(
            "get_active_nodes", +[](const expression<T> &instance) { return v_to_l(instance.get_active_nodes()); },
            "Gets the idx of the active nodes in the current chromosome")
        .def("get_n", &expression<T>::get_n, "Gets the number of inputs of the dCGP expression")
        .def("get_m", &expression<T>::get_m, "Gets the number of outputs of the dCGP expression")
        .def("get_rows", &expression<T>::get_r, "Gets the number of rows of the dCGP expression")
        .def("get_cols", &expression<T>::get_c, "Gets the number of columns of the dCGP expression")
        .def("get_levels_back", &expression<T>::get_l, "Gets the number of levels-back allowed for the dCGP expression")
        .def(
            "get_arity", +[](const expression<T> &instance) { return v_to_l(instance.get_arity()); },
            "get_arity()\nget_arity(node_id)\nGets the arity of the basis functions of the dCGP expression. Either "
            "the whole vector or that of a single node.")
        .def(
            "get_arity", +[](const expression<T> &instance, unsigned node_id) { return instance.get_arity(node_id); },
            (bp::arg("node_id")))
        .def(
            "get_gene_idx", +[](const expression<T> &instance) { return v_to_l(instance.get_gene_idx()); },
            "get_gene_idx()\nGets a vector containing the indexes in the chromosome where each node starts to be "
            "expressed.")
        .def(
            "get_f", +[](const expression<T> &instance) { return v_to_l(instance.get_f()); },
            "Gets the kernel functions")
        .def(
            "mutate", +[](expression<T> &instance, const bp::object &in) { instance.mutate(l_to_v<unsigned>(in)); },
            expression_mutate_doc().c_str(), bp::arg("idxs"))
        .def("mutate_random", &expression<T>::mutate_random,
             "mutate_random(N = 1)\nMutates N randomly selected genes within its allowed bounds", bp::arg("N"))
        .def("mutate_active", &expression<T>::mutate_active,
             "mutate_active(N = 1)\nMutates N randomly selected active genes within their allowed bounds",
             (bp::arg("N") = 1))
        .def("mutate_active_cgene", &expression<T>::mutate_active_cgene,
             "mutate_active_cgene(N = 1)\nMutates N randomly selected active connections within their allowed bounds",
             (bp::arg("N") = 1))
        .def("mutate_ogene", &expression<T>::mutate_ogene,
             "mutate_ogene(N = 1)\nMutates N randomly selected output genes connection within their allowed bounds",
             (bp::arg("N") = 1))
        .def(
            "mutate_active_fgene", &expression<T>::mutate_active_fgene,
            "mutate_active_fgene(N = 1)\nMutates N randomly selected active function genes within their allowed bounds",
            (bp::arg("N") = 1))
        // The parallelism for the loss computation is switched off in python as pitonic kernels can produce a crash
        .def(
            "loss",
            +[](const expression<T> &instance, const bp::object &points, const bp::object &labels,
                const std::string &loss) {
                auto parallel = 0u;
                return instance.loss(to_vv<T>(points), to_vv<T>(labels), loss, parallel);
            },
            expression_loss_doc().c_str(), (bp::arg("points"), bp::arg("labels"), bp::arg("loss")))
        .add_property(
            "eph_val", +[](const expression<T> &instance) { return v_to_l(instance.get_eph_val()); },
            +[](expression<T> &instance, const bp::object &eph_val) { instance.set_eph_val(l_to_v<T>(eph_val)); })
        .add_property(
            "eph_symb", +[](const expression<T> &instance) { return v_to_l(instance.get_eph_symb()); },
            +[](expression<T> &instance, const bp::object &eph_symb) {
                instance.set_eph_symb(l_to_v<std::string>(eph_symb));
            });
}

template <typename T>
void expose_expression_weighted(std::string type)
{
    std::string class_name = "expression_weighted_" + type;
    bp::class_<expression_weighted<T>, bp::bases<expression<T>>>(class_name.c_str(), bp::no_init)
        // Constructor with seed
        .def("__init__",
             bp::make_constructor(
                 +[](unsigned in, unsigned out, unsigned rows, unsigned cols, unsigned levelsback,
                     const bp::object &arity, const bp::object &kernels, unsigned seed) {
                     auto kernels_v = l_to_v<kernel<T>>(kernels);
                     bp::extract<unsigned> is_int(arity);
                     if (is_int.check()) { // arity is passed as an integer
                         unsigned ar = bp::extract<unsigned>(arity);
                         return ::new expression_weighted<T>(in, out, rows, cols, levelsback, ar, kernels_v, seed);
                     } else { // arity is passed as something else, a list is assumed
                         auto varity = l_to_v<unsigned>(arity);
                         return ::new expression_weighted<T>(in, out, rows, cols, levelsback, varity, kernels_v, seed);
                     }
                 },
                 bp::default_call_policies(),
                 (bp::arg("inputs"), bp::arg("outputs"), bp::arg("rows"), bp::arg("cols"), bp::arg("levels_back"),
                  bp::arg("arity"), bp::arg("kernels"), bp::arg("seed"))),
             expression_init_doc(type).c_str())
        // Constructor with no seed
        .def("__init__",
             bp::make_constructor(
                 +[](unsigned in, unsigned out, unsigned rows, unsigned cols, unsigned levelsback,
                     const bp::object &arity, const bp::object &kernels) {
                     auto kernels_v = l_to_v<kernel<T>>(kernels);
                     bp::extract<unsigned> is_int(arity);
                     if (is_int.check()) { // arity is passed as an integer
                         unsigned ar = bp::extract<unsigned>(arity);
                         return ::new expression_weighted<T>(in, out, rows, cols, levelsback, ar, kernels_v,
                                                             std::random_device()());
                     } else { // arity is passed as something else, a list is assumed
                         auto varity = l_to_v<unsigned>(arity);
                         return ::new expression_weighted<T>(in, out, rows, cols, levelsback, varity, kernels_v,
                                                             std::random_device()());
                     }
                 },
                 bp::default_call_policies(),
                 (bp::arg("inputs"), bp::arg("outputs"), bp::arg("rows"), bp::arg("cols"), bp::arg("levels_back"),
                  bp::arg("arity"), bp::arg("kernels"))),
             expression_init_doc(type).c_str())
        .def(
            "__repr__",
            +[](const expression_weighted<T> &instance) -> std::string {
                std::ostringstream oss;
                oss << instance;
                return oss.str();
            })
        .def(
            "__call__",
            +[](const expression_weighted<T> &instance, const bp::object &in) {
                try {
                    auto v = l_to_v<T>(in);
                    return v_to_l(instance(v));
                } catch (...) {
                    PyErr_Clear();
                    auto v = l_to_v<std::string>(in);
                    return v_to_l(instance(v));
                }
            })
        .def("set_weight", &expression_weighted<T>::set_weight, expression_weighted_set_weight_doc().c_str(),
             (bp::arg("node_id"), bp::arg("input_id"), bp::arg("weight")))
        .def(
            "set_weights",
            +[](expression_weighted<T> &instance, const bp::object &weights) {
                instance.set_weights(l_to_v<T>(weights));
            },
            expression_weighted_set_weights_doc().c_str(), (bp::arg("weights")))
        .def("get_weight", &expression_weighted<T>::get_weight, expression_weighted_get_weight_doc().c_str(),
             (bp::arg("node_id"), bp::arg("input_id")))
        .def(
            "get_weights", +[](expression_weighted<T> &instance) { return v_to_l(instance.get_weights()); },
            "Gets all weights");
}

template <typename T>
void expose_expression_ann(std::string type)
{
    std::string class_name = "expression_ann_" + type;
    bp::class_<expression_ann, bp::bases<expression<T>>>(class_name.c_str(), bp::no_init)
        // Constructor with seed
        .def("__init__",
             bp::make_constructor(
                 +[](unsigned in, unsigned out, unsigned rows, unsigned cols, unsigned levelsback,
                     const bp::object &arity, const bp::object &kernels, unsigned seed) {
                     auto kernels_v = l_to_v<kernel<T>>(kernels);
                     bp::extract<unsigned> is_int(arity);
                     if (is_int.check()) { // arity is passed as an integer
                         unsigned ar = bp::extract<unsigned>(arity);
                         return ::new expression_ann(in, out, rows, cols, levelsback, ar, kernels_v, seed);
                     } else { // arity is passed as something else, a list is assumed
                         auto varity = l_to_v<unsigned>(arity);
                         return ::new expression_ann(in, out, rows, cols, levelsback, varity, kernels_v, seed);
                     }
                 },
                 bp::default_call_policies(),
                 (bp::arg("inputs"), bp::arg("outputs"), bp::arg("rows"), bp::arg("cols"), bp::arg("levels_back"),
                  bp::arg("arity"), bp::arg("kernels"), bp::arg("seed"))),
             expression_init_doc(type).c_str())
        // Constructor with no seed
        .def("__init__",
             bp::make_constructor(
                 +[](unsigned in, unsigned out, unsigned rows, unsigned cols, unsigned levelsback,
                     const bp::object &arity, const bp::object &kernels) {
                     auto kernels_v = l_to_v<kernel<double>>(kernels);
                     bp::extract<unsigned> is_int(arity);
                     if (is_int.check()) { // arity is passed as an integer
                         unsigned ar = bp::extract<unsigned>(arity);
                         return ::new expression_ann(in, out, rows, cols, levelsback, ar, kernels_v,
                                                     std::random_device()());
                     } else { // arity is passed as something else, a list is assumed
                         auto varity = l_to_v<unsigned>(arity);
                         return ::new expression_ann(in, out, rows, cols, levelsback, varity, kernels_v,
                                                     std::random_device()());
                     }
                 },
                 bp::default_call_policies(),
                 (bp::arg("inputs"), bp::arg("outputs"), bp::arg("rows"), bp::arg("cols"), bp::arg("levels_back"),
                  bp::arg("arity"), bp::arg("kernels"))),
             expression_init_doc(type).c_str())
        .def(
            "__repr__",
            +[](const expression_ann &instance) -> std::string {
                std::ostringstream oss;
                oss << instance;
                return oss.str();
            })
        .def(
            "__call__",
            +[](const expression_ann &instance, const bp::object &in) {
                try {
                    auto v = l_to_v<double>(in);
                    return v_to_l(instance(v));
                } catch (...) {
                    PyErr_Clear();
                    auto v = l_to_v<std::string>(in);
                    return v_to_l(instance(v));
                }
            })
        .def("set_bias", &expression_ann::set_bias, expression_ann_set_bias_doc().c_str(),
             (bp::arg("node_id"), bp::arg("bias")))
        .def(
            "set_biases",
            +[](expression_ann &instance, const bp::object &biases) { instance.set_biases(l_to_v<double>(biases)); },
            expression_ann_set_biases_doc().c_str(), (bp::arg("biases")))
        .def("get_bias", &expression_ann::get_bias, expression_ann_get_bias_doc().c_str(), (bp::arg("node_id")))
        .def(
            "get_biases", +[](expression_ann &instance) { return v_to_l(instance.get_biases()); }, "Gets all biases")
        .def(
            "set_weight", +[](expression_ann &instance, unsigned idx, double w) { instance.set_weight(idx, w); },
            (bp::arg("idx"), bp::arg("value")))
        .def(
            "set_weight",
            +[](expression_ann &instance, unsigned node_id, unsigned input_id, double w) {
                instance.set_weight(node_id, input_id, w);
            },
            expression_ann_set_weight_doc().c_str(), (bp::arg("node_id"), bp::arg("input_id"), bp::arg("value")))
        .def(
            "set_weights",
            +[](expression_ann &instance, const bp::object &weights) { instance.set_weights(l_to_v<T>(weights)); },
            expression_weighted_set_weights_doc().c_str(), (bp::arg("weights")))
        .def("set_output_f", &expression_ann::set_output_f, expression_ann_set_output_f_doc().c_str(),
             (bp::arg("f_id")))
        .def(
            "get_weight", +[](expression_ann &instance, unsigned idx) { return instance.get_weight(idx); },
            (bp::arg("idx")))
        .def(
            "get_weight",
            +[](expression_ann &instance, unsigned node_id, unsigned input_id) {
                return instance.get_weight(node_id, input_id);
            },
            expression_ann_get_weight_doc().c_str(), (bp::arg("node_id"), bp::arg("input_id")))
        .def(
            "get_weights", +[](expression_ann &instance) { return v_to_l(instance.get_weights()); }, "Gets all weights")
        .def("n_active_weights", &expression_ann::n_active_weights, expression_ann_n_active_weights_doc().c_str(),
             bp::arg("unique") = false)
        .def(
            "randomise_weights",
            +[](expression_ann &instance, double mean, double std, unsigned seed) {
                return instance.randomise_weights(mean, std, seed);
            },
            expression_ann_randomise_weights_doc().c_str(),
            (bp::arg("mean") = 0., bp::arg("std") = 0.1, bp::arg("seed")))
        .def(
            "randomise_weights",
            +[](expression_ann &instance, double mean, double std) {
                return instance.randomise_weights(mean, std, std::random_device()());
            },
            (bp::arg("mean") = 0., bp::arg("std") = 0.1))
        .def(
            "randomise_biases",
            +[](expression_ann &instance, double mean, double std, unsigned seed) {
                return instance.randomise_biases(mean, std, seed);
            },
            expression_ann_randomise_biases_doc().c_str(),
            (bp::arg("mean") = 0., bp::arg("std") = 0.1, bp::arg("seed")))
        .def(
            "randomise_biases",
            +[](expression_ann &instance, double mean, double std) {
                return instance.randomise_biases(mean, std, std::random_device()());
            },
            (bp::arg("mean") = 0., bp::arg("std") = 0.1))
        .def(
            "sgd",
            +[](expression_ann &instance, const bp::object &points, const bp::object &labels, double l_rate,
                unsigned batch_size, const std::string &loss, unsigned parallel, bool shuffle) {
                auto d = to_vv<double>(points);
                auto l = to_vv<double>(labels);
                return instance.sgd(d, l, l_rate, batch_size, loss, parallel, shuffle);
            },
            expression_ann_sgd_doc().c_str(),
            (bp::arg("points"), bp::arg("labels"), bp::arg("lr"), bp::arg("batch_size"), bp::arg("loss"),
             bp::arg("parallel") = 0u, bp::arg("shuffle") = true));
}

using gym_ptr = void (*)(std::vector<std::vector<double>> &, std::vector<std::vector<double>> &);
template <gym_ptr F>
inline void expose_data_from_the_gym(const std::string &name, const std::string &docstring = std::string{"ds"})
{
    bp::def(
        name.c_str(),
        +[]() {
            std::vector<std::vector<double>> points, labels;
            F(points, labels);
            return bp::make_tuple(vvector_to_ndarr<double>(points), vvector_to_ndarr<double>(labels));
        },
        docstring.c_str());
}

BOOST_PYTHON_MODULE(core)
{
    bp::docstring_options doc_options;

    // Init numpy.
    // NOTE: only the second import is strictly necessary. We run a first import from BP
    // because that is the easiest way to detect whether numpy is installed or not (rather
    // than trying to figure out a way to detect it from import_array()).
    try {
        bp::import("numpy.core.multiarray");
    } catch (...) {
        dcgpy::builtin().attr("print")(
            u8"\033[91m====ERROR====\nThe NumPy module could not be imported. "
            u8"Please make sure that NumPy has been correctly installed.\n====ERROR====\033[0m");
        throw std::runtime_error("Import error");
    }
    dcgpy::numpy_import_array();

    doc_options.enable_all();
    doc_options.disable_cpp_signatures();
    doc_options.disable_py_signatures();

    expose_kernel<double>("double");
    expose_kernel_set<double>("double");
    expose_expression<double>("double");
    expose_expression_weighted<double>("double");
    expose_expression_ann<double>("double");

    expose_kernel<gdual_d>("gdual_double");
    expose_kernel_set<gdual_d>("gdual_double");
    expose_expression<gdual_d>("gdual_double");
    expose_expression_weighted<gdual_d>("gdual_double");

    expose_kernel<gdual_v>("gdual_vdouble");
    expose_kernel_set<gdual_v>("gdual_vdouble");
    expose_expression<gdual_v>("gdual_vdouble");
    expose_expression_weighted<gdual_v>("gdual_vdouble");

    // Making data from the gym available in python
    expose_data_from_the_gym<&gym::generate_koza_quintic>("generate_koza_quintic", generate_koza_quintic_doc());
    // From Our paper
    expose_data_from_the_gym<&gym::generate_P1>("generate_P1", generate_P1_doc());
    expose_data_from_the_gym<&gym::generate_P2>("generate_P2", generate_P2_doc());
    expose_data_from_the_gym<&gym::generate_P3>("generate_P3", generate_P3_doc());
    expose_data_from_the_gym<&gym::generate_P4>("generate_P4", generate_P4_doc());
    expose_data_from_the_gym<&gym::generate_P5>("generate_P5", generate_P5_doc());
    expose_data_from_the_gym<&gym::generate_P6>("generate_P6", generate_P6_doc());
    expose_data_from_the_gym<&gym::generate_P7>("generate_P7", generate_P7_doc());
    // From Vladi paper
    expose_data_from_the_gym<&gym::generate_kotanchek>("generate_kotanchek", generate_kotanchek_doc());
    expose_data_from_the_gym<&gym::generate_salutowicz>("generate_salutowicz", generate_salutowicz_doc());
    expose_data_from_the_gym<&gym::generate_salutowicz2d>("generate_salutowicz2d", generate_salutowicz2d_doc());
    expose_data_from_the_gym<&gym::generate_uball5d>("generate_uball5d", generate_uball5d_doc());
    expose_data_from_the_gym<&gym::generate_ratpol3d>("generate_ratpol3d", generate_ratpol3d_doc());
    expose_data_from_the_gym<&gym::generate_sinecosine>("generate_sinecosine", generate_sinecosine_doc());
    expose_data_from_the_gym<&gym::generate_vladi7>("generate_vladi7");
    expose_data_from_the_gym<&gym::generate_vladi8>("generate_vladi8");
    // NIST data
    expose_data_from_the_gym<&gym::generate_chwirut1>("generate_chwirut1");
    expose_data_from_the_gym<&gym::generate_chwirut2>("generate_chwirut2");
    expose_data_from_the_gym<&gym::generate_daniel_wood>("generate_daniel_wood");
    expose_data_from_the_gym<&gym::generate_gauss1>("generate_gauss1");
    expose_data_from_the_gym<&gym::generate_kirby2>("generate_kirby2");
    expose_data_from_the_gym<&gym::generate_lanczos2>("generate_lanczos2");
    expose_data_from_the_gym<&gym::generate_misra1b>("generate_misra1b");

    // Define a cleanup functor to be run when the module is unloaded.
    struct dcgp_cleanup_functor {
        void operator()() const
        {
            std::cout << "Shutting down the thread pool.\n";
            piranha::thread_pool_shutdown<void>();
        }
    };
    // Expose it.
    bp::class_<dcgp_cleanup_functor> cl_c("_dcgp_cleanup_functor", bp::init<>());
    cl_c.def("__call__", &dcgp_cleanup_functor::operator());
    // Register it.
    bp::object atexit_mod = bp::import("atexit");
    atexit_mod.attr("register")(dcgp_cleanup_functor{});
}
