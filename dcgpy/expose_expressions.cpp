#include <boost/numeric/conversion/cast.hpp>
#include <pybind11/pybind11.h>
#include <sstream>
#include <string>
#include <vector>

#include <dcgp/expression.hpp>
#include <dcgp/expression_ann.hpp>
#include <dcgp/expression_weighted.hpp>
#include <dcgp/kernel.hpp>

#include "docstrings.hpp"

using namespace dcgp;
using namespace audi;
namespace py = pybind11;

namespace dcgpy
{

template <typename T>
void expose_expression(const py::module &m, std::string type)
{
    std::string class_name = "expression_" + type;
    auto exp_ = py::class_<expression<T>>(m, class_name.c_str(), "A CGP expression");
    exp_.def(py::init<>())
        .def(py::init<unsigned, unsigned, unsigned, unsigned, unsigned, std::vector<unsigned>, std::vector<kernel<T>>,
                      unsigned, unsigned>(),
             py::arg("inputs"), py::arg("outputs"), py::arg("rows"), py::arg("cols"), py::arg("levels_back"),
             py::arg("arity"), py::arg("kernels"), py::arg("n_eph"), py::arg("seed"), expression_init_doc(type).c_str())
        .def(py::init<unsigned, unsigned, unsigned, unsigned, unsigned, std::vector<unsigned>, std::vector<kernel<T>>,
                      unsigned>(),
             py::arg("inputs"), py::arg("outputs"), py::arg("rows"), py::arg("cols"), py::arg("levels_back"),
             py::arg("arity"), py::arg("kernels"), py::arg("n_eph"), expression_init_doc(type).c_str())
        .def(
            "__repr__",
            +[](const expression<T> &instance) -> std::string {
                std::ostringstream oss;
                oss << instance;
                return oss.str();
            })
        .def(
            "__call__", [](const expression<T> &instance, const std::vector<T> &v) { return instance(v); },
            "Call operator from values")
        .def(
            "__call__", [](const expression<T> &instance, const std::vector<std::string> &v) { return instance(v); },
            "Call operator from strings")
        .def("set", &expression<T>::set, expression_set_doc().c_str(), py::arg("chromosome"))
        .def("set_f_gene", &expression<T>::set_f_gene, expression_set_f_gene_doc().c_str(), py::arg("node_id"),
             py::arg("f_id"))
        .def("get", &expression<T>::get, "Gets the expression chromosome")
        .def("get_lb", &expression<T>::get_lb, "Gets the lower bounds of the chromosome")
        .def("get_ub", &expression<T>::get_ub, "Gets the upper bounds of the chromosome")
        .def("get_active_genes", &expression<T>::get_active_genes,
             "Gets the idx of the active genes in the current chromosome (numbering is from 0)")
        .def("get_active_nodes", &expression<T>::get_active_nodes,
             "Gets the idx of the active nodes in the current chromosome")
        .def("get_n", &expression<T>::get_n, "Gets the number of inputs of the dCGP expression")
        .def("get_m", &expression<T>::get_m, "Gets the number of outputs of the dCGP expression")
        .def("get_rows", &expression<T>::get_r, "Gets the number of rows of the dCGP expression")
        .def("get_cols", &expression<T>::get_c, "Gets the number of columns of the dCGP expression")
        .def("get_levels_back", &expression<T>::get_l, "Gets the number of levels-back allowed for the dCGP expression")
        .def(
            "get_arity", [](const expression<T> &instance) { return instance.get_arity(); },
            "get_arity()\nget_arity(node_id)\nGets the arity of the basis functions of the dCGP expression. Either "
            "the whole vector or that of a single node.")
        .def(
            "get_arity", [](const expression<T> &instance, unsigned node_id) { return instance.get_arity(node_id); },
            py::arg("node_id"))
        .def("get_gene_idx", &expression<T>::get_gene_idx,
             "get_gene_idx()\nGets a vector containing the indexes in the chromosome where each node starts to be "
             "expressed.")
        .def("get_f", &expression<T>::get_f, "Gets the kernel functions")
        .def(
            "mutate", [](expression<T> &in, unsigned idx) { return in.mutate(idx); }, expression_mutate_doc().c_str(),
            py::arg("idxs"))
        .def(
            "mutate", [](expression<T> &in, const std::vector<unsigned> &idxs) { return in.mutate(idxs); },
            expression_mutate_doc().c_str(), py::arg("idxs"))
        .def("mutate_random", &expression<T>::mutate_random,
             "mutate_random(N = 1)\nMutates N randomly selected genes within its allowed bounds", py::arg("N"))
        .def("mutate_active", &expression<T>::mutate_active,
             "mutate_active(N = 1)\nMutates N randomly selected active genes within their allowed bounds",
             py::arg("N") = 1)
        .def("mutate_active_cgene", &expression<T>::mutate_active_cgene,
             "mutate_active_cgene(N = 1)\nMutates N randomly selected active connections within their allowed bounds",
             py::arg("N") = 1)
        .def("mutate_ogene", &expression<T>::mutate_ogene,
             "mutate_ogene(N = 1)\nMutates N randomly selected output genes connection within their allowed bounds",
             py::arg("N") = 1)
        .def(
            "mutate_active_fgene", &expression<T>::mutate_active_fgene,
            "mutate_active_fgene(N = 1)\nMutates N randomly selected active function genes within their allowed bounds",
            py::arg("N") = 1)
        // The parallelism for the loss computation is switched off in python as pythonic kernels would
        // produce a crash if evaluated in multiple threads.
        .def(
            "loss",
            [](const expression<T> &instance, const std::vector<std::vector<T>> &points,
               const std::vector<std::vector<T>> &labels,
               const std::string &loss) { return instance.loss(points, labels, loss, 0u); },
            expression_loss_doc().c_str(), py::arg("points"), py::arg("labels"), py::arg("loss"))
        .def_property("eph_val", &expression<T>::get_eph_val, &expression<T>::set_eph_val)
        .def_property("eph_symb", &expression<T>::get_eph_symb, &expression<T>::set_eph_symb);
}

// template <typename T>
// void expose_expression_weighted(std::string type)
//{
//    std::string class_name = "expression_weighted_" + type;
//    bp::class_<expression_weighted<T>, bp::bases<expression<T>>>(class_name.c_str(), bp::no_init)
//        // Constructor with seed
//        .def("__init__",
//             bp::make_constructor(
//                 +[](unsigned in, unsigned out, unsigned rows, unsigned cols, unsigned levelsback,
//                     const bp::object &arity, const bp::object &kernels, unsigned seed) {
//                     auto kernels_v = l_to_v<kernel<T>>(kernels);
//                     bp::extract<unsigned> is_int(arity);
//                     if (is_int.check()) { // arity is passed as an integer
//                         unsigned ar = bp::extract<unsigned>(arity);
//                         return ::new expression_weighted<T>(in, out, rows, cols, levelsback, ar, kernels_v, seed);
//                     } else { // arity is passed as something else, a list is assumed
//                         auto varity = l_to_v<unsigned>(arity);
//                         return ::new expression_weighted<T>(in, out, rows, cols, levelsback, varity, kernels_v,
//                         seed);
//                     }
//                 },
//                 bp::default_call_policies(),
//                 (bp::arg("inputs"), bp::arg("outputs"), bp::arg("rows"), bp::arg("cols"), bp::arg("levels_back"),
//                  bp::arg("arity"), bp::arg("kernels"), bp::arg("seed"))),
//             expression_init_doc(type).c_str())
//        // Constructor with no seed
//        .def("__init__",
//             bp::make_constructor(
//                 +[](unsigned in, unsigned out, unsigned rows, unsigned cols, unsigned levelsback,
//                     const bp::object &arity, const bp::object &kernels) {
//                     auto kernels_v = l_to_v<kernel<T>>(kernels);
//                     bp::extract<unsigned> is_int(arity);
//                     if (is_int.check()) { // arity is passed as an integer
//                         unsigned ar = bp::extract<unsigned>(arity);
//                         return ::new expression_weighted<T>(in, out, rows, cols, levelsback, ar, kernels_v,
//                                                             std::random_device()());
//                     } else { // arity is passed as something else, a list is assumed
//                         auto varity = l_to_v<unsigned>(arity);
//                         return ::new expression_weighted<T>(in, out, rows, cols, levelsback, varity, kernels_v,
//                                                             std::random_device()());
//                     }
//                 },
//                 bp::default_call_policies(),
//                 (bp::arg("inputs"), bp::arg("outputs"), bp::arg("rows"), bp::arg("cols"), bp::arg("levels_back"),
//                  bp::arg("arity"), bp::arg("kernels"))),
//             expression_init_doc(type).c_str())
//        .def(
//            "__repr__",
//            +[](const expression_weighted<T> &instance) -> std::string {
//                std::ostringstream oss;
//                oss << instance;
//                return oss.str();
//            })
//        .def(
//            "__call__",
//            +[](const expression_weighted<T> &instance, const bp::object &in) {
//                try {
//                    auto v = l_to_v<T>(in);
//                    return v_to_l(instance(v));
//                } catch (...) {
//                    PyErr_Clear();
//                    auto v = l_to_v<std::string>(in);
//                    return v_to_l(instance(v));
//                }
//            })
//        .def("set_weight", &expression_weighted<T>::set_weight, expression_weighted_set_weight_doc().c_str(),
//             (bp::arg("node_id"), bp::arg("input_id"), bp::arg("weight")))
//        .def(
//            "set_weights",
//            +[](expression_weighted<T> &instance, const bp::object &weights) {
//                instance.set_weights(l_to_v<T>(weights));
//            },
//            expression_weighted_set_weights_doc().c_str(), (bp::arg("weights")))
//        .def("get_weight", &expression_weighted<T>::get_weight, expression_weighted_get_weight_doc().c_str(),
//             (bp::arg("node_id"), bp::arg("input_id")))
//        .def(
//            "get_weights", +[](expression_weighted<T> &instance) { return v_to_l(instance.get_weights()); },
//            "Gets all weights");
//}
//
// template <typename T>
// void expose_expression_ann(std::string type)
//{
//    std::string class_name = "expression_ann_" + type;
//    bp::class_<expression_ann, bp::bases<expression<T>>>(class_name.c_str(), bp::no_init)
//        // Constructor with seed
//        .def("__init__",
//             bp::make_constructor(
//                 +[](unsigned in, unsigned out, unsigned rows, unsigned cols, unsigned levelsback,
//                     const bp::object &arity, const bp::object &kernels, unsigned seed) {
//                     auto kernels_v = l_to_v<kernel<T>>(kernels);
//                     bp::extract<unsigned> is_int(arity);
//                     if (is_int.check()) { // arity is passed as an integer
//                         unsigned ar = bp::extract<unsigned>(arity);
//                         return ::new expression_ann(in, out, rows, cols, levelsback, ar, kernels_v, seed);
//                     } else { // arity is passed as something else, a list is assumed
//                         auto varity = l_to_v<unsigned>(arity);
//                         return ::new expression_ann(in, out, rows, cols, levelsback, varity, kernels_v, seed);
//                     }
//                 },
//                 bp::default_call_policies(),
//                 (bp::arg("inputs"), bp::arg("outputs"), bp::arg("rows"), bp::arg("cols"), bp::arg("levels_back"),
//                  bp::arg("arity"), bp::arg("kernels"), bp::arg("seed"))),
//             expression_init_doc(type).c_str())
//        // Constructor with no seed
//        .def("__init__",
//             bp::make_constructor(
//                 +[](unsigned in, unsigned out, unsigned rows, unsigned cols, unsigned levelsback,
//                     const bp::object &arity, const bp::object &kernels) {
//                     auto kernels_v = l_to_v<kernel<double>>(kernels);
//                     bp::extract<unsigned> is_int(arity);
//                     if (is_int.check()) { // arity is passed as an integer
//                         unsigned ar = bp::extract<unsigned>(arity);
//                         return ::new expression_ann(in, out, rows, cols, levelsback, ar, kernels_v,
//                                                     std::random_device()());
//                     } else { // arity is passed as something else, a list is assumed
//                         auto varity = l_to_v<unsigned>(arity);
//                         return ::new expression_ann(in, out, rows, cols, levelsback, varity, kernels_v,
//                                                     std::random_device()());
//                     }
//                 },
//                 bp::default_call_policies(),
//                 (bp::arg("inputs"), bp::arg("outputs"), bp::arg("rows"), bp::arg("cols"), bp::arg("levels_back"),
//                  bp::arg("arity"), bp::arg("kernels"))),
//             expression_init_doc(type).c_str())
//        .def(
//            "__repr__",
//            +[](const expression_ann &instance) -> std::string {
//                std::ostringstream oss;
//                oss << instance;
//                return oss.str();
//            })
//        .def(
//            "__call__",
//            +[](const expression_ann &instance, const bp::object &in) {
//                try {
//                    auto v = l_to_v<double>(in);
//                    return v_to_l(instance(v));
//                } catch (...) {
//                    PyErr_Clear();
//                    auto v = l_to_v<std::string>(in);
//                    return v_to_l(instance(v));
//                }
//            })
//        .def("set_bias", &expression_ann::set_bias, expression_ann_set_bias_doc().c_str(),
//             (bp::arg("node_id"), bp::arg("bias")))
//        .def(
//            "set_biases",
//            +[](expression_ann &instance, const bp::object &biases) { instance.set_biases(l_to_v<double>(biases)); },
//            expression_ann_set_biases_doc().c_str(), (bp::arg("biases")))
//        .def("get_bias", &expression_ann::get_bias, expression_ann_get_bias_doc().c_str(), (bp::arg("node_id")))
//        .def(
//            "get_biases", +[](expression_ann &instance) { return v_to_l(instance.get_biases()); }, "Gets all biases")
//        .def(
//            "set_weight", +[](expression_ann &instance, unsigned idx, double w) { instance.set_weight(idx, w); },
//            (bp::arg("idx"), bp::arg("value")))
//        .def(
//            "set_weight",
//            +[](expression_ann &instance, unsigned node_id, unsigned input_id, double w) {
//                instance.set_weight(node_id, input_id, w);
//            },
//            expression_ann_set_weight_doc().c_str(), (bp::arg("node_id"), bp::arg("input_id"), bp::arg("value")))
//        .def(
//            "set_weights",
//            +[](expression_ann &instance, const bp::object &weights) { instance.set_weights(l_to_v<T>(weights)); },
//            expression_weighted_set_weights_doc().c_str(), (bp::arg("weights")))
//        .def("set_output_f", &expression_ann::set_output_f, expression_ann_set_output_f_doc().c_str(),
//             (bp::arg("f_id")))
//        .def(
//            "get_weight", +[](expression_ann &instance, unsigned idx) { return instance.get_weight(idx); },
//            (bp::arg("idx")))
//        .def(
//            "get_weight",
//            +[](expression_ann &instance, unsigned node_id, unsigned input_id) {
//                return instance.get_weight(node_id, input_id);
//            },
//            expression_ann_get_weight_doc().c_str(), (bp::arg("node_id"), bp::arg("input_id")))
//        .def(
//            "get_weights", +[](expression_ann &instance) { return v_to_l(instance.get_weights()); }, "Gets all
//            weights")
//        .def("n_active_weights", &expression_ann::n_active_weights, expression_ann_n_active_weights_doc().c_str(),
//             bp::arg("unique") = false)
//        .def(
//            "randomise_weights",
//            +[](expression_ann &instance, double mean, double std, unsigned seed) {
//                return instance.randomise_weights(mean, std, seed);
//            },
//            expression_ann_randomise_weights_doc().c_str(),
//            (bp::arg("mean") = 0., bp::arg("std") = 0.1, bp::arg("seed")))
//        .def(
//            "randomise_weights",
//            +[](expression_ann &instance, double mean, double std) {
//                return instance.randomise_weights(mean, std, std::random_device()());
//            },
//            (bp::arg("mean") = 0., bp::arg("std") = 0.1))
//        .def(
//            "randomise_biases",
//            +[](expression_ann &instance, double mean, double std, unsigned seed) {
//                return instance.randomise_biases(mean, std, seed);
//            },
//            expression_ann_randomise_biases_doc().c_str(),
//            (bp::arg("mean") = 0., bp::arg("std") = 0.1, bp::arg("seed")))
//        .def(
//            "randomise_biases",
//            +[](expression_ann &instance, double mean, double std) {
//                return instance.randomise_biases(mean, std, std::random_device()());
//            },
//            (bp::arg("mean") = 0., bp::arg("std") = 0.1))
//        .def(
//            "sgd",
//            +[](expression_ann &instance, const bp::object &points, const bp::object &labels, double l_rate,
//                unsigned batch_size, const std::string &loss, unsigned parallel, bool shuffle) {
//                auto d = to_vv<double>(points);
//                auto l = to_vv<double>(labels);
//                return instance.sgd(d, l, l_rate, batch_size, loss, parallel, shuffle);
//            },
//            expression_ann_sgd_doc().c_str(),
//            (bp::arg("points"), bp::arg("labels"), bp::arg("lr"), bp::arg("batch_size"), bp::arg("loss"),
//             bp::arg("parallel") = 0u, bp::arg("shuffle") = true));

void expose_expressions(const py::module &m)
{
    // double
    expose_expression<double>(m, "double");
    // expose_expression_weighted<double>(m, "double");
    // expose_expression_ann<double>(m, "double");
    // gdual_d
    expose_expression<gdual_d>(m, "gdual_double");
    // expose_expression_weighted<gdual_d>("gdual_double");
    // gdual_v
    expose_expression<gdual_v>(m, "gdual_vdouble");
    // expose_expression_weighted<gdual_v>("gdual_vdouble");
}
} // namespace dcgpy