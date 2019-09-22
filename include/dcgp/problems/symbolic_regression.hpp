#ifndef DCGP_SYMBOLIC_REGRESSION_H
#define DCGP_SYMBOLIC_REGRESSION_H
#include <pagmo/types.hpp>
#include <vector>

namespace dcgp
{
class symbolic_regression
{
public:
    // Constructor from std::vectors
    symbolic_regression(const std::vector<std::vector<double>> &points, const std::vector<std::vector<double>> &labels)
    {
        // We check that points is not empty
        if (points.size() == 0) {
            throw std::invalid_argument("The size of the input data (points) is zero.");
        }
        // We check that labels are defined for each input point
        if (points.size() != labels.size()) {
            throw std::invalid_argument("The number of input data (points) is " + std::to_string(points.size())
                                        + " while the number of labels is " + std::to_string(labels.size())
                                        + ". They should be equal.");
        }
        auto p_size = points[0].size();
        if (!std::all_of(points.begin(), points.end(),
                         [p_size](const std::vector<double> &p) { return p.size() == p_size; })) {
            throw std::invalid_argument("The input data (points) is inconsistent: all points must have the same "
                                        "dimension, while I detect differences.");
        }
        if (!std::all_of(labels.begin(), labels.end(),
                         [p_size](const std::vector<double> &p) { return p.size() == p_size; })) {
            throw std::invalid_argument("The labels are inconsistent: all labels must have the same "
                                        "dimension, while I detect differences.");
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
