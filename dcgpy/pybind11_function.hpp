#ifndef DCGPY_FUNCTION_INNER_FOR_PYOBJECT_HPP
#define DCGPY_FUNCTION_INNER_FOR_PYOBJECT_HPP

#include <audi/audi.hpp>
#include <pybind11/pybind11.h>
#include <sstream>
#include <string>
#include <vector>

#include <dcgp/s11n.hpp>

#include "common_utils.hpp"
#include "docstrings.hpp"

using namespace dcgp;
using namespace dcgpy;
using namespace audi;
namespace py = pybind11;

namespace dcgp::detail
{

template <typename R, typename... Args>
struct function_inner<py::object, R, Args...> final : function_inner_base<R, Args...> {
    // We just need the def ctor, delete everything else.
    function_inner() = default;
    function_inner(const function_inner &) = delete;
    function_inner(function_inner &&) = delete;
    function_inner &operator=(const function_inner &) = delete;
    function_inner &operator=(function_inner &&) = delete;

    // Constructor from generic python object.
    explicit function_inner(const py::object &o)
    {
        m_value = dcgpy::deepcopy(o);
    }

    // The clone method, used in the copy constructor of function.
    virtual std::unique_ptr<function_inner_base<R, Args...>> clone() const override final
    {
        return std::make_unique<function_inner>(m_value);
    }

    // Mandatory methods.
    virtual R operator()(Args &&... args) const override final
    {
        return  py::cast<R>(m_value(std::forward<Args>(args)...));
    }

    virtual pagmo::thread_safety get_thread_safety() const override final
    {
        return pagmo::thread_safety::none;
    }

    template <typename Archive>
    void save(Archive &ar, unsigned) const
    {
        ar << boost::serialization::base_object<function_inner_base<R, Args...>>(*this);
        ar << dcgpy::object_to_vchar(m_value);
    }
    template <typename Archive>
    void load(Archive &ar, unsigned)
    {
        ar >> boost::serialization::base_object<function_inner_base<R, Args...>>(*this);
        std::vector<char> v;
        ar >> v;
        m_value = dcgpy::vchar_to_object(v);
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()

    py::object m_value;
};

} // namespace dcgp::detail

#endif