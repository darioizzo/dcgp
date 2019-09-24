#ifndef DCGP_SYMBOLIC_REGRESSION_H
#define DCGP_SYMBOLIC_REGRESSION_H
#include <boost/numeric/conversion/cast.hpp>
#include <boost/range/algorithm/transform.hpp>
#include <pagmo/io.hpp>
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
    /// Default constructor
    /**
     * A default constructor is needed by the pagmo UDP interface, but it should not be used.
     * It constructs a list of 1 empty points/labels vector and a dummy cgp member.
     * It is then guaranteed that m_points[0] and m_labels[0] can be accessed.
     */
    symbolic_regression() : m_points(1), m_labels(1), m_cgp(1u, 1u, 1u, 1u, 1u, 2u, kernel_set<double>({"sum"})(), 0u)
    {
    }

    /// Constructor
    /**
     * Constructs a symbolic_regression optimization problem compatible with the pagmo UDP interface.
     *
     * @param[in] points number of inputs (independent variables).
     * @param[in] labels number of outputs (dependent variables).
     */
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
        // We initialize the dcgp expression
        kernel_set<double> basic_set({"sum", "diff", "mul", "div"});
        m_cgp = expression<double>(n, m, 2, 10, 10, 2, basic_set(), random_device::next());
    }
    /// Fitness computation
    /**
     * Computes the fitness for this UDP
     *
     * @param x the decision vector.
     *
     * @return the fitness of \p x.
     */
    pagmo::vector_double fitness(const pagmo::vector_double &x) const
    {
        // We need to make a copy of the chromosome as to represents its genes as unsigned
        std::vector<unsigned> xu(x.size());
        std::transform(x.begin(), x.end(), xu.begin(), [](double a) { return boost::numeric_cast<unsigned>(a); });
        m_cgp.set(xu);
        // We initialize the fitness
        std::vector<double> f(1, 0);
        // We compute the MSE loss splitting the data in n batches (if possible).
        f[0] = m_cgp.loss(m_points, m_labels, "MSE", 0u);
        return f;
    }
    /// Box-bounds
    /**
     *
     * It returns the box-bounds for this UDP.
     *
     * @return the lower and upper bounds for each of the decision vector components
     */
    std::pair<pagmo::vector_double, pagmo::vector_double> get_bounds() const
    {
        std::vector<double> lb(m_cgp.get_lb().begin(), m_cgp.get_lb().end());
        std::vector<double> ub(m_cgp.get_ub().begin(), m_cgp.get_ub().end());
        return {lb, ub};
    }
    /// Integer dimension
    /**
     * It returns the integer dimension of the problem.
     *
     * @return the integer dimension of the problem.
     */
    pagmo::vector_double::size_type get_nix() const
    {
        return m_cgp.get_lb().size();
    }
    /// Problem name
    /**
     * @return a string containing the problem name
     */
    std::string get_name() const
    {
        return "CGP symbolic regressor";
    }
    /// Extra info
    /**
     * @return a string containing extra problem information.
     */
    std::string get_extra_info() const
    {
        std::ostringstream ss;
        pagmo::stream(ss, "\tInput dimension: ", m_points[0].size(), "\n");
        pagmo::stream(ss, "\tOutput dimension: ", m_labels[0].size(), "\n");
        pagmo::stream(ss, "\tData size: ", m_points.size(), "\n");
        pagmo::stream(ss, "\tKernels: ", m_cgp.get_f(), "\n");
        return ss.str();
    }

    /// Human-readable representation of a decision vector.
    /**
     * @param[in] x a valid chromosome.
     *
     * @return a string containing the mathematical expression represented by *x*.
     */
    std::string pretty(const pagmo::vector_double &x) const
    {
        // We need to make a copy of the chromosome as to represents its genes as unsigned
        std::vector<unsigned> xu(x.size());
        std::transform(x.begin(), x.end(), xu.begin(), [](double a) { return boost::numeric_cast<unsigned>(a); });
        m_cgp.set(xu);
        std::ostringstream ss;
        std::vector<std::string> symbols;
        for (decltype(m_points[0].size()) i = 0u; i < m_points[0].size(); ++i) {
            symbols.push_back("x" + std::to_string(i));
        }
        pagmo::stream(ss, m_cgp(symbols));
        return ss.str();
    }

private:
    std::vector<std::vector<double>> m_points;
    std::vector<std::vector<double>> m_labels;
    mutable expression<double> m_cgp;
};
} // namespace dcgp
#endif
