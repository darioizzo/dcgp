#ifndef DCGP_DETAIL_S11N_WRAPPERS_HPP
#define DCGP_DETAIL_S11N_WRAPPERS_HPP

#include <utility>

namespace dcgp
{

namespace detail
{

// A few helpers to give a cereal-like syntax
// to Boost.serialization.
template <typename Archive>
inline void archive(Archive &)
{
}

template <typename Archive, typename Arg0, typename... Args>
inline void archive(Archive &ar, Arg0 &&arg0, Args &&... args)
{
    ar &std::forward<Arg0>(arg0);
    archive(ar, std::forward<Args>(args)...);
}

template <typename Archive>
inline void to_archive(Archive &)
{
}

template <typename Archive, typename Arg0, typename... Args>
inline void to_archive(Archive &ar, Arg0 &&arg0, Args &&... args)
{
    ar << std::forward<Arg0>(arg0);
    to_archive(ar, std::forward<Args>(args)...);
}

template <typename Archive>
inline void from_archive(Archive &)
{
}

template <typename Archive, typename Arg0, typename... Args>
inline void from_archive(Archive &ar, Arg0 &&arg0, Args &&... args)
{
    ar >> std::forward<Arg0>(arg0);
    from_archive(ar, std::forward<Args>(args)...);
}

} // namespace detail

} // namespace dcgp

#endif
