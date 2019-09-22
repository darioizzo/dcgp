#ifndef DCGP_SYMBOLIC_REGRESSION_H
#define DCGP_SYMBOLIC_REGRESSION_H
#include <pagmo/types.hpp>
#include <vector>

#include <dcgp/expression.hpp>
#include <dcgp/kernel_set.hpp>
#include <dcgp/rng.hpp>

namespace dcgp
{
class symbolic_regression
{
public:
    // Constructor from std::vectors
    explicit symbolic_regression(const std::vector<std::vector<double>> &points,
                                 const std::vector<std::vector<double>> &labels)
        : m_points(points), m_labels(labels), m_cgp(1u, 1u, 1u, 1u, 1u, 2u, kernel_set<double>({"sum"})(), 0u)
    {
        // 1 - We check that points is not an empty vector.
        if (points.size() == 0) {
            throw std::invalid_argument("The size of the input data (points) is zero.");
        }
        // 2 - We check labels and points have the same (non-empty) size
        if (points.size() != labels.size()) {
            throw std::invalid_argument("The number of input data (points) is " + std::to_string(points.size())
                                        + " while the number of labels is " + std::to_string(labels.size())
                                        + ". They should be equal.");
        }
        // 3 - We check that all p in points have the same size
        unsigned n = static_cast<unsigned>(points[0].size());
        if (!std::all_of(points.begin(), points.end(), [n](const std::vector<double> &p) { return p.size() == n; })) {
            throw std::invalid_argument("The input data (points) is inconsistent: all points must have the same "
                                        "dimension, while I detect differences.");
        }
        // 4 - We check that all l in labels have the same size
        unsigned m = static_cast<unsigned>(labels[0].size());
        if (!std::all_of(labels.begin(), labels.end(), [m](const std::vector<double> &l) { return l.size() == m; })) {
            throw std::invalid_argument("The labels are inconsistent: all labels must have the same "
                                        "dimension, while I detect differences.");
        }
        kernel_set<double> basic_set({"sum", "diff", "mul", "div"});
        m_cgp = expression<double>(n, m, 2, 10, 10, 2, basic_set(), random_device::next());
    }
    // Specified by the pagmo UDP interface.
    pagmo::vector_double fitness(const pagmo::vector_double & x) const
    {
        std::vector<unsigned> xu(x.begin(), x.end());
        m_cgp.set(xu);
        std::vector<double> f(1,0);
        f[0] = m_cgp.loss(m_points, m_labels, "MSE", true);
        return f;
    }
    // Specified by the pagmo UDP interface.
    std::pair<pagmo::vector_double, pagmo::vector_double> get_bounds() const
    {
        std::vector<double> lb(m_cgp.get_lb().begin(), m_cgp.get_lb().end());
        std::vector<double> ub(m_cgp.get_ub().begin(), m_cgp.get_ub().end());
        return {lb, ub};
    }

private:
    std::vector<std::vector<double>> m_points;
    std::vector<std::vector<double>> m_labels;
    mutable expression<double> m_cgp;
};
} // namespace dcgp
#endif
