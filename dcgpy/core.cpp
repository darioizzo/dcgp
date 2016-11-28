#include <audi/audi.hpp>
#include <boost/python.hpp>
#include <vector>
#include <string>
#include <sstream>
#include <functional> //std::function

#include "common_utils.hpp"
#include "../include/kernel.hpp"
#include "../include/kernel_set.hpp"
#include "../include/expression.hpp"
#include "../include/expression_weighted.hpp"
#include "docstrings.hpp"

using namespace dcgp;
using namespace dcgpy;
using namespace audi;
namespace bp = boost::python;

template <typename T>
void expose_kernel(const std::string &type)
{
    std::string class_name = "kernel_" + type;
    bp::class_<kernel<T>>(class_name.c_str(), bp::no_init)
    .def("__init__", bp::make_constructor(
        +[](const bp::object &obj1, const bp::object &obj2, const std::string &name)
            {
                std::function<T(const std::vector<T>&)> my_function = [obj1](const std::vector<T>& x)
                {
                    T in = bp::extract<T>(obj1(v_to_l(x)));
                    return in;
                };
                std::function<std::string(const std::vector<std::string>&)> my_print_function = [obj2](const std::vector<std::string>& x)
                {
                    std::string in = bp::extract<std::string>(obj2(v_to_l(x)));
                    return in;
                };
                return ::new kernel<T>(my_function, my_print_function, name);
            },
            bp::default_call_policies(),
            (bp::arg("callable_f"), bp::arg("callable_s"), bp::arg("name"))
            ),
            kernel_init_doc(type).c_str()
        )
    .def("__call__",
        +[](kernel<T> &instance, const bp::object &in)
        {
            try {
                auto v = l_to_v<T>(in);
                return bp::object(instance(v));
            } catch (...) {
                PyErr_Clear();
                auto v = l_to_v<std::string>(in);
                return bp::object(instance(v));
            }
        }
    )
    .def("__repr__",
        +[](const kernel<T> &instance) -> std::string
        {
            std::ostringstream oss;
            oss << instance;
            return oss.str();
        }
    );
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
    bp::class_<kernel_set<T>>(class_name.c_str(), bp::no_init)
    .def("__init__", bp::make_constructor(
        +[](const bp::object &obj1)
        {
            auto a = l_to_v<std::string>(obj1);
            return ::new kernel_set<T>(a);
        },
        bp::default_call_policies(),
        (bp::arg("kernels"))
        ),
        kernel_set_init_doc(type).c_str()
    )
    .def("__call__",
        +[](kernel_set<T> &instance)
        {
            return v_to_l(instance());
        }
    )
    .def("__repr__",
        +[](const kernel_set<T> &instance) -> std::string
        {
            std::ostringstream oss;
            oss << instance;
            return oss.str();
        }
    )
    .def("push_back", (void (kernel_set<T>::*)(std::string)) &kernel_set<T>::push_back, "Adds one more kernel to the set by common name")
    .def("push_back", (void (kernel_set<T>::*)(const kernel<T>&)) &kernel_set<T>::push_back, "Adds one more kernel to the set")
    .def( "__getitem__", &wrap_operator<T>);
}

template <typename T>
void expose_expression(std::string type)
{
    std::string class_name = "expression_" + type;
    bp::class_<expression<T>>(class_name.c_str(), bp::no_init)
    .def("__init__", bp::make_constructor(
        +[](unsigned int in, unsigned int out, unsigned int rows, unsigned int cols, unsigned int levelsback, unsigned int arity, const bp::object &kernels, unsigned int seed)
        {
            auto kernels_v = l_to_v<kernel<T>>(kernels);
            return ::new expression<T>(in, out, rows, cols, levelsback, arity, kernels_v, seed);
        },
        bp::default_call_policies(),
        (bp::arg("inputs"),bp::arg("outputs"),bp::arg("rows"),bp::arg("cols"),bp::arg("levels_back"),bp::arg("arity"),bp::arg("kernels"),bp::arg("seed"))
        ),
        expression_init_doc(type).c_str()
    )
    .def("__repr__",
        +[](const expression<T> &instance) -> std::string
        {
            std::ostringstream oss;
            oss << instance;
            return oss.str();
        }
    )
    .def("__call__",
        +[](const expression<T> &instance, const bp::object &in)
        {
            try {
                auto v = l_to_v<T>(in);
                return v_to_l(instance(v));
            } catch (...) {
                PyErr_Clear();
                auto v = l_to_v<std::string>(in);
                return v_to_l(instance(v));
            }
        }
    )
    .def("set", +[](expression<T> &instance, const bp::object &in){instance.set(l_to_v<unsigned int>(in));}, "Sets the expression chromosme", bp::arg("chromosme"))
    .def("get", +[](const expression<T> &instance){return v_to_l(instance.get());}, "Gets the expression chromosme")
    .def("get_lb", +[](const expression<T> &instance){return v_to_l(instance.get_lb());}, "Gets the lower bounds of the chromosome")
    .def("get_ub", +[](const expression<T> &instance){return v_to_l(instance.get_ub());}, "Gets the upper bounds of the chromosome")
    .def("get_active_genes", +[](const expression<T> &instance){return v_to_l(instance.get_active_genes());}, "Gets the idx of the active genes in the current chromosome (numbering is from 0)")
    .def("get_active_nodes", +[](const expression<T> &instance){return v_to_l(instance.get_active_nodes());}, "Gets the idx of the active nodes in the current chromosome.")
    .def("get_n", &expression<T>::get_n, "Gets the number of inputs of the d_CGP expression")
    .def("get_m", &expression<T>::get_m, "Gets the number of outputs of the d_CGP expression")
    .def("get_rows", &expression<T>::get_rows, "Gets the number of rows of the d_CGP expression")
    .def("get_cols", &expression<T>::get_cols, "Gets the number of columns of the d_CGP expression")
    .def("get_levels_back", &expression<T>::get_levels_back, "Gets the number of levels-back allowed for the d_CGP expression")
    .def("get_arity", &expression<T>::get_arity, "Gets the arity of the basis functions of the d_CGP expression")
    .def("get_f", +[](const expression<T> &instance){return v_to_l(instance.get_f());}, "Gets the kernel functions")
    .def("mutate", +[](expression<T> &instance, const bp::object &in){instance.mutate(l_to_v<unsigned int>(in));}, "Mutates the selected genes within its allowed bounds.", bp::arg("idxs"))
    .def("mutate_random", &expression<T>::mutate_random, "Mutates N randomly selected genes within its allowed bounds.")
    .def("mutate_active", &expression<T>::mutate_active, "Mutates N randomly selected active genes within its allowed bounds.")
    .def("mutate_active_cgene", &expression<T>::mutate_active_cgene, "Mutates exactly one randomly selected active connection genes within its allowed bounds.")
    .def("mutate_ogene", &expression<T>::mutate_ogene, "Mutates exactly one randomly selected output genes within its allowed bounds.")
    .def("mutate_active_fgene", &expression<T>::mutate_active_fgene, "Mutates exactly one randomly selected active function genes within its allowed bounds.");
}

template <typename T>
void expose_expression_weighted(std::string type)
{
    std::string class_name = "expression_weighted_" + type;
    bp::class_<expression_weighted<T>, bp::bases<expression<T> > >(class_name.c_str(), bp::no_init)
    .def("__init__", bp::make_constructor(
        +[](unsigned int in, unsigned int out, unsigned int rows, unsigned int cols, unsigned int levelsback, unsigned int arity, const bp::object &kernels, unsigned int seed)
        {
            auto kernels_v = l_to_v<kernel<T>>(kernels);
            return ::new expression_weighted<T>(in, out, rows, cols, levelsback, arity, kernels_v, seed);
        },
        bp::default_call_policies(),
        (bp::arg("in"),bp::arg("out"),bp::arg("rows"),bp::arg("cols"),bp::arg("levels_back"),bp::arg("arity"),bp::arg("kernels"),bp::arg("seed"))
        ),
        expression_init_doc(type).c_str()
    )
    .def("__repr__",
        +[](const expression_weighted<T> &instance) -> std::string
        {
            std::ostringstream oss;
            oss << instance;
            return oss.str();
        }
    )
    .def("__call__",
        +[](const expression_weighted<T> &instance, const bp::object &in)
        {
            try {
                auto v = l_to_v<T>(in);
                return v_to_l(instance(v));
            } catch (...) {
                PyErr_Clear();
                auto v = l_to_v<std::string>(in);
                return v_to_l(instance(v));
            }
        }
    )
    .def("set_weight", &expression_weighted<T>::set_weight, "Sets a weight", (bp::arg("node_id"), bp::arg("input_id"), bp::arg("weight")))
    .def("set_weights",
        +[] (expression_weighted<T> &instance, const bp::object &weights)
        {
            instance.set_weights(l_to_v<T>(weights));
        },
        "Sets all weights at once",
        (bp::arg("weights"))
    )
    .def("get_weights",
        +[] (expression_weighted<T> &instance)
        {
            return v_to_l(instance.get_weights());
        },
        "Gets all weights"
    );
}

BOOST_PYTHON_MODULE(_core)
{
    expose_kernel<double>("double");
    expose_kernel_set<double>("double");
    expose_expression<double>("double");
    expose_expression_weighted<double>("double");
    expose_kernel<gdual_d>("gdual_double");
    expose_kernel_set<gdual_d>("gdual_double");
    expose_expression<gdual_d>("gdual_double");
    expose_expression_weighted<gdual_d>("gdual_double");
    expose_kernel<gdual_v>("gdual_vdouble");
    expose_kernel_set<gdual_v>("gdual_vdouble");
    expose_expression<gdual_v>("gdual_vdouble");
    expose_expression_weighted<gdual_v>("gdual_vdouble");

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
