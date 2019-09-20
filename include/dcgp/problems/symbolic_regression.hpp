#ifndef DCGP_SYMBOLIC_REGRESSION_H
#include <pagmo/types.hpp>

namespace dcgp
{
class symbolic_regression
{
public:
    pagmo::vector_double fitness(const pagmo::vector_double &) const
    {
        return {1.};
    }
    std::pair<pagmo::vector_double, pagmo::vector_double> get_bounds() const
    {
        return {{0., 1.}, {2., 3.}};
    }
};
} // namespace dcgp
#endif
