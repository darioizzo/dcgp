#ifndef DCGP_FUNCTION_H
#define DCGP_FUNCTION_H

#include <memory>
#include <stdexcept>
#include <type_traits>
#include <utility>

#include <boost/type_traits/integral_constant.hpp>
#include <boost/type_traits/is_virtual_base_of.hpp>

#include <pagmo/threading.hpp>
#include <pagmo/type_traits.hpp>

#include <dcgp/s11n.hpp>

// NOTE: we disable address tracking for all user-defined classes. The reason is that even if the final
// classes (i.e., function) use value semantics, the internal implementation details use old-style
// OO construct (i.e., base classes, pointers, etc.). By default, Boost serialization wants to track
// the addresses of these internal implementation-detail classes, and this has some undesirable consequences
// (for instance, when deserializing a function object in a variable and then moving it into another
// one, which is a pattern we sometimes use in order to provide increased exception safety).
//
// See also:
// https://www.boost.org/doc/libs/1_70_0/libs/serialization/doc/special.html#objecttracking
// https://www.boost.org/doc/libs/1_70_0/libs/serialization/doc/traits.html#level
#define DCGP_S11N_FUNCTION_EXPORT_KEY(id, func, ...)                                                                   \
    namespace dcgp::s11n_names                                                                                         \
    {                                                                                                                  \
    using udf_##id = dcgp::detail::function_inner<func, __VA_ARGS__>;                                                  \
    }                                                                                                                  \
    BOOST_CLASS_EXPORT_KEY2(dcgp::s11n_names::udf_##id, "udf " #id)                                                    \
    BOOST_CLASS_TRACKING(dcgp::s11n_names::udf_##id, boost::serialization::track_never)

#define DCGP_S11N_FUNCTION_IMPLEMENT(id, func, ...)                                                                    \
    namespace dcgp::s11n_names                                                                                         \
    {                                                                                                                  \
    using udf_##id = dcgp::detail::function_inner<func, __VA_ARGS__>;                                                  \
    }                                                                                                                  \
    BOOST_CLASS_EXPORT_IMPLEMENT(dcgp::s11n_names::udf_##id)

#define DCGP_S11N_FUNCTION_EXPORT(id, func, ...)                                                                       \
    DCGP_S11N_FUNCTION_EXPORT_KEY(id, func, __VA_ARGS__)                                                               \
    DCGP_S11N_FUNCTION_IMPLEMENT(id, func, __VA_ARGS__)

namespace dcgp::detail
{

template <typename R, typename... Args>
struct function_inner_base {
    virtual ~function_inner_base() {}
    virtual std::unique_ptr<function_inner_base> clone() const = 0;
    virtual R operator()(Args &&...) const = 0;
    virtual pagmo::thread_safety get_thread_safety() const = 0;
    template <typename Archive>
    void serialize(Archive &, unsigned)
    {
    }
};

template <typename T, typename R, typename... Args>
struct function_inner final : function_inner_base<R, Args...> {
    // We just need the def ctor, delete everything else.
    function_inner() = default;
    function_inner(const function_inner &) = delete;
    function_inner(function_inner &&) = delete;
    function_inner &operator=(const function_inner &) = delete;
    function_inner &operator=(function_inner &&) = delete;

    // Constructors from T (copy and move variants).
    explicit function_inner(const T &x) : m_value(x) {}
    explicit function_inner(T &&x) : m_value(std::move(x)) {}

    // The clone method, used in the copy constructor of function.
    virtual std::unique_ptr<function_inner_base<R, Args...>> clone() const override final
    {
        return std::make_unique<function_inner>(m_value);
    }

    // Mandatory methods.
    virtual R operator()(Args &&... args) const override final
    {
        return m_value(std::forward<Args>(args)...);
    }

    // Optional methods.
    virtual pagmo::thread_safety get_thread_safety() const override final
    {
        return get_thread_safety_impl(m_value);
    }
    template <typename U, std::enable_if_t<pagmo::has_get_thread_safety<U>::value, int> = 0>
    static pagmo::thread_safety get_thread_safety_impl(const U &value)
    {
        return value.get_thread_safety();
    }
    template <typename U, std::enable_if_t<!pagmo::has_get_thread_safety<U>::value, int> = 0>
    static pagmo::thread_safety get_thread_safety_impl(const U &)
    {
        return pagmo::thread_safety::basic;
    }

    // Serialization.
    template <typename Archive>
    void serialize(Archive &ar, unsigned)
    {
        ar &boost::serialization::base_object<function_inner_base<R, Args...>>(*this);
        ar &m_value;
    }

    T m_value;
};

} // namespace dcgp::detail

namespace boost
{

template <typename T, typename R, typename... Args>
struct is_virtual_base_of<dcgp::detail::function_inner_base<R, Args...>, dcgp::detail::function_inner<T, R, Args...>>
    : false_type {
};

} // namespace boost

namespace dcgp
{

template <typename R, typename... Args>
class function;

template <typename R, typename... Args>
class function<R(Args...)>
{
    // Helpful alias.
    template <typename T>
    using uncvref_t = std::remove_cv_t<std::remove_reference_t<T>>;

    // Function pointer type returning R and taking Args
    // as parameters.
    using f_ptr = R (*)(Args...);

    // Dispatching for the generic ctor. We have a special case if T is
    // a function type, in which case we will manually do the conversion to
    // function pointer and delegate to the other overload.
    template <typename T>
    explicit function(T &&x, std::true_type) : function(static_cast<f_ptr>(std::forward<T>(x)), std::false_type{})
    {
    }
    template <typename T>
    explicit function(T &&x, std::false_type)
        : m_ptr(std::make_unique<detail::function_inner<uncvref_t<T>, R, Args...>>(std::forward<T>(x)))
    {
    }

public:
    function() : function(static_cast<f_ptr>(nullptr)) {}
    function(const function &other) : m_ptr(other.ptr()->clone()), m_thread_safety(other.m_thread_safety) {}
    function(function &&) noexcept = default;
    function &operator=(function &&other) noexcept
    {
        if (this != &other) {
            m_ptr = std::move(other.m_ptr);
            m_thread_safety = std::move(other.m_thread_safety);
        }
        return *this;
    }
    function &operator=(const function &other)
    {
        // Copy ctor + move assignment.
        return *this = function(other);
    }

    template <typename T, std::enable_if_t<!std::is_same_v<function, uncvref_t<T>>, int> = 0>
    function(T &&x) : function(std::forward<T>(x), std::is_function<uncvref_t<T>>{})
    {
        // Assign the thread safety level.
        m_thread_safety = ptr()->get_thread_safety();
    }

    // Extraction and related.
    template <typename T>
    const T *extract() const noexcept
    {
        auto p = dynamic_cast<const detail::function_inner<T, R, Args...> *>(ptr());
        return p == nullptr ? nullptr : &(p->m_value);
    }
    template <typename T>
    T *extract() noexcept
    {
        auto p = dynamic_cast<detail::function_inner<T, R, Args...> *>(ptr());
        return p == nullptr ? nullptr : &(p->m_value);
    }
    template <typename T>
    bool is() const noexcept
    {
        return extract<T>() != nullptr;
    }

    R operator()(Args... args) const
    {
        auto ptr = extract<f_ptr>();
        if (ptr != nullptr && *ptr == nullptr) {
            throw std::runtime_error(
                "This dcgp::function object cannot be invoked because it contains a null pointer to a C++ function");
        }
        return m_ptr->operator()(std::forward<Args>(args)...);
    }

    // Thread safety level.
    pagmo::thread_safety get_thread_safety() const
    {
        return m_thread_safety;
    }

    // Check if the function was not
    // moved-from.
    bool is_valid() const
    {
        return static_cast<bool>(m_ptr);
    }

    // Serialisation support.
    template <typename Archive>
    void save(Archive &ar, unsigned) const
    {
        ar << m_ptr;
        ar << m_thread_safety;
    }
    template <typename Archive>
    void load(Archive &ar, unsigned)
    {
        // Deserialize in a separate object and move it in later, for exception safety.
        function tmp_func;
        ar >> tmp_func.m_ptr;
        ar >> tmp_func.m_thread_safety;
        *this = std::move(tmp_func);
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()

private:
    // Just two small helpers to make sure that whenever we require
    // access to the pointer it actually points to something.
    detail::function_inner_base<R, Args...> const *ptr() const
    {
        assert(m_ptr.get() != nullptr);
        return m_ptr.get();
    }
    detail::function_inner_base<R, Args...> *ptr()
    {
        assert(m_ptr.get() != nullptr);
        return m_ptr.get();
    }

private:
    // Pointer to the inner base function.
    std::unique_ptr<detail::function_inner_base<R, Args...>> m_ptr;
    // Thread safety.
    pagmo::thread_safety m_thread_safety;
};

} // namespace dcgp

// Disable class tracking for all instances of dcgp::function.
namespace boost::serialization
{

template <typename R, typename... Args>
struct tracking_level<dcgp::function<R, Args...>> {
    typedef mpl::integral_c_tag tag;
    typedef mpl::int_<track_never> type;
    BOOST_STATIC_CONSTANT(int, value = tracking_level::type::value);
    BOOST_STATIC_ASSERT(
        (mpl::greater<implementation_level<dcgp::function<R, Args...>>, mpl::int_<primitive_type>>::value));
};

} // namespace boost::serialization

#endif
