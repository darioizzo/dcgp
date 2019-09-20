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
        if (point.size() != this->get_n()) {
            throw std::invalid_argument("When computing the loss the point dimension (input) seemed wrong, it was: "
                                        + std::to_string(point.size())
                                        + " while I expected: " + std::to_string(this->get_n()));
        }
        if (prediction.size() != this->get_m()) {
            throw std::invalid_argument(
                "When computing the loss the prediction dimension (output) seemed wrong, it was: "
                + std::to_string(prediction.size()) + " while I expected: " + std::to_string(this->get_m()));
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
