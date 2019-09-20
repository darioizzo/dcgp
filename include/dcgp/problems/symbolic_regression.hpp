#ifndef DCGP_SYMBOLIC_REGRESSION_H
#include <pagmo/types.hpp>
#include <vector.hpp>

namespace dcgp
{
class symbolic_regression
{
public:
    // Constructor from std::vectors
    symbolic_regression(const std::vector<std::vector<double>> &points, const std::vector<std::vector<double>> &labels)
    {
        if (points.size() != labels.size()) {
            throw std::invalid_argument("The number of input data (points) is " + std::to_string(point.size())
                                        + " while the number of labels is " + std::to_string(labels.size())
                                        + ". They should be equal.");
        }
        if !std::all_of(points.begin(), points.end(), [](const std::vector<double> &p) { return p.size == p_size; }) {

        } 
    }
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
