#ifndef DCGP_SYMBOLIC_REGRESSION_H
#define DCGP_SYMBOLIC_REGRESSION_H
#include <audi/gdual.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/range/algorithm/transform.hpp>
#include <pagmo/io.hpp>
#include <pagmo/population.hpp>
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
     * It is guaranteed that m_points[0] and m_labels[0] can be accessed.
     */
    symbolic_regression()
        : m_points(1), m_labels(1), m_r(1), m_c(1), m_l(1), m_arity(2), m_f(kernel_set<double>({"sum"})()),
          m_parallel_batches(0u)
    {
    }

    /// Constructor
    /**
     * Constructs a symbolic_regression optimization problem compatible with the pagmo UDP interface.
     *
     * @param[in] points number of inputs (independent variables).
     * @param[in] labels number of outputs (dependent variables).
     * @param[in] r number of rows of the dCGP.
     * @param[in] c number of columns of the dCGP.
     * @param[in] l number of levels-back allowed in the dCGP.
     * @param[in] arity arity of the basis functions.
     * @param[in] f function set. An std::vector of dcgp::kernel<expression::type>.
     * @param[in] parallel_batches number of parallel batches.
     *
     * @throws std::invalid_argument if points and labels are not consistent.
     * @throws std::invalid_argument if the CGP related parameters (i.e. *r*, *c*, etc...) are malformed.
     */
    symbolic_regression(const std::vector<std::vector<double>> &points, const std::vector<std::vector<double>> &labels,
                        unsigned r = 1,     // n. rows
                        unsigned c = 10,    // n. columns
                        unsigned l = 11,    // n. levels-back
                        unsigned arity = 2, // basis functions' arity
                        std::vector<kernel<double>> f
                        = kernel_set<double>({"sum", "diff", "mul", "pdiv"})(), // functions
                        unsigned n_eph = 0u,                                    // number of ephemeral constants
                        unsigned parallel_batches = 0u                          // number of parallel batches
                        )
        : m_points(points), m_labels(labels), m_r(r), m_c(c), m_l(l), m_arity(arity), m_f(f), m_n_eph(n_eph),
          m_parallel_batches(parallel_batches)
    {
        unsigned n;
        unsigned m;
        // We check the inputs.
        sanity_checks(n, m);
        // We initialize the cgp expression
        auto seed = random_device::next();
        m_cgp = expression<double>(n, m, m_r, m_c, m_l, m_arity, m_f, m_n_eph, seed);
        // We initialize the dcgp expression
        kernel_set<audi::gdual_d> f_g;
        for (const auto &ker : f) {
            auto name = ker.get_name();
            // If protected division is used in the cgp the unprotected one will be used in the dCGP
            // as to keep derivative sanity
            if (name.compare("pdiv") == 0) {
                f_g.push_back("div");
            } else {
                f_g.push_back(ker.get_name());
            }
        }
        m_dcgp = expression<gdual_d>(n, m, m_r, m_c, m_l, m_arity, f_g(), m_n_eph, seed);
        // We initialize the dpoints/dduals
        m_dpoints.clear();
        m_dlabels.clear();
        for (const auto &point : m_points) {
            std::vector<audi::gdual_d> point_gdual;
            for (const auto &item : point) {
                point_gdual.emplace_back(item);
            }
            m_dpoints.push_back(point_gdual);
        }
        for (const auto &label : m_labels) {
            std::vector<audi::gdual_d> label_gdual;
            for (const auto &item : label) {
                label_gdual.emplace_back(item);
            }
            m_dlabels.push_back(label_gdual);
        }
        // We create the symbol set of the differentials here for efficiency.
        // They are used in the gradient computation.
        for (const auto &symb : m_dcgp.get_eph_symb()) {
            m_deph_symb.push_back("d" + symb);
        }
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
        if (x == m_cache.first) {
            return m_cache.second;
        } else {
            std::vector<double> retval(1, 0);
            // The chromosome has a floating point part (the ephemeral constants) and an integer part (the encoded CGP).
            // 1 - We extract the integer part and represent it as an unsigned vector to set the CGP expression.
            std::vector<unsigned> xu(x.size() - m_n_eph);
            std::transform(x.data() + m_n_eph, x.data() + x.size(), xu.begin(),
                           [](double a) { return boost::numeric_cast<unsigned>(a); });
            m_cgp.set(xu);
            // 2 - We set the floating point part as ephemeral constants.
            std::vector<double> eph_val(x.data(), x.data() + m_n_eph);
            m_cgp.set_eph_val(eph_val);
            // 3 - We compute the MSE loss splitting the data in n batches.
            // TODO: make this work also when m_parallel_batches does not divide exactly the data size.
            retval[0] = m_cgp.loss(m_points, m_labels, "MSE", m_parallel_batches);
            return retval;
        }
    }

    /// Gradient computation
    /**
     * Computes the gradient for the continuous part of this UDP
     *
     * @param x the decision vector.
     *
     * @return the gradient of \p x.
     */
    pagmo::vector_double gradient(const pagmo::vector_double &x) const
    {
        std::vector<double> retval(m_n_eph, 0);
        // The chromosome has a floating point part (the ephemeral constants) and an integer part (the encoded CGP).
        // 1 - We extract the integer part and represent it as an unsigned vector to set the CGP expression.
        std::vector<unsigned> xu(x.size() - m_n_eph);
        std::transform(x.data() + m_n_eph, x.data() + x.size(), xu.begin(),
                       [](double a) { return boost::numeric_cast<unsigned>(a); });
        m_dcgp.set(xu);
        // 2 - We use the floating point part of the chromosome to set ephemeral constants.
        std::vector<audi::gdual_d> eph_val;
        for (decltype(m_n_eph) i = 0u; i < m_n_eph; ++i) {
            eph_val.emplace_back(x[i], m_dcgp.get_eph_symb()[i], 1u); // Only first derivative is needed
        }
        m_dcgp.set_eph_val(eph_val);
        // 3 - We compute the MSE loss splitting the data in n batches.
        // TODO: make this work also when m_parallel_batches does not divide exactly the data size.
        auto loss = m_dcgp.loss(m_dpoints, m_dlabels, "MSE", m_parallel_batches);
        // since we have also computed the loss value, we store it in a cache so that fitness can be called and
        // not cause reavaluation of a cgp.
        m_cache.first = x;
        m_cache.second.resize(1);
        m_cache.second[0] = loss.constant_cf();
        loss.extend_symbol_set(m_deph_symb);
        if (!(loss.get_order() == 0u)) { // this happens when input terminals of the eph constants are inactive
                                         // (gradient is then zero)
            for (decltype(m_n_eph) i = 0u; i < m_n_eph; ++i) {
                std::vector<double> coeff(m_n_eph, 0.);
                coeff[i] = 1.;
                retval[i] = loss.get_derivative(coeff);
            }
        }
        return retval;
    }

    /// Sparsity pattern
    /**
     * Returns the sparsity pattern. The sparsity patter is dense in the continuous part of the chromosome.
     * We assume all eph constants are in the expression. If not, we deal with a vanishing gradient later.
     *
     * @return the sparsity pattern.
     */
    pagmo::sparsity_pattern gradient_sparsity() const
    {
        pagmo::sparsity_pattern retval;
        for (decltype(m_n_eph) i = 0u; i < m_n_eph; ++i) {
            retval.push_back({0u, i});
        }
        return retval;
    }

    /// Box-bounds
    /**
     *
     * Returns the box-bounds for this UDP.
     *
     * @return the lower and upper bounds for each of the decision vector components
     */
    std::pair<pagmo::vector_double, pagmo::vector_double> get_bounds() const
    {
        // Bounds on ephemeral constants are -10, 10;
        std::vector<double> lb(m_cgp.get_lb().size() + m_n_eph, -10.);
        std::vector<double> ub(m_cgp.get_ub().size() + m_n_eph, 10.);
        // Bounds on the CGP encoding are derived from the dcgp::expression
        std::copy(m_cgp.get_lb().begin(), m_cgp.get_lb().end(), lb.data() + m_n_eph);
        std::copy(m_cgp.get_ub().begin(), m_cgp.get_ub().end(), ub.data() + m_n_eph);
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
        return m_cgp.get_lb().size() - m_n_eph;
    }

    /// Problem name
    /**
     * @return a string containing the problem name
     */
    std::string get_name() const
    {
        return "a CGP symbolic regression problem";
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
        std::vector<unsigned> xu(x.size() - m_n_eph);
        std::transform(x.data() + m_n_eph, x.data() + x.size(), xu.data(),
                       [](double a) { return boost::numeric_cast<unsigned>(a); });
        m_cgp.set(xu);
        // 2 - We set the floating point part as ephemeral constants.
        std::vector<double> eph_val(x.data(), x.data() + m_n_eph);
        m_cgp.set_eph_val(eph_val);

        std::ostringstream ss;
        std::vector<std::string> symbols;
        for (decltype(m_points[0].size()) i = 0u; i < m_points[0].size(); ++i) {
            symbols.push_back("x" + std::to_string(i));
        }
        pagmo::stream(ss, m_cgp(symbols));
        return ss.str();
    }

    /// Gets the inner CGP
    /**
     * The access to the inner CGP is offered in the public interface to allow evolve methods in UDAs
     * to reuse the same object and perform mutations via it. This is a hack to interface pagmo with
     * dCGP. Alternatives would be a friendship relation (uughhh) or construct a new CGP object within
     * the evolve each time (seems expensive). So here it is, FOR USE ONLY IN udas::evolve methods.
     */
    const expression<double> &get_cgp() const
    {
        return m_cgp;
    }

private:
    inline void sanity_checks(unsigned &n, unsigned &m) const
    {
        // 1 - We check that points is not an empty vector.
        if (m_points.size() == 0) {
            throw std::invalid_argument("The size of the input data (points) is zero.");
        }
        n = static_cast<unsigned>(m_points[0].size());
        m = static_cast<unsigned>(m_labels[0].size());
        // 2 - We check labels and points have the same (non-empty) size
        if (m_points.size() != m_labels.size()) {
            throw std::invalid_argument("The number of input data (points) is " + std::to_string(m_points.size())
                                        + " while the number of labels is " + std::to_string(m_labels.size())
                                        + ". They should be equal.");
        }
        // 3 - We check that all p in points have the same size
        if (!std::all_of(m_points.begin(), m_points.end(),
                         [n](const std::vector<double> &p) { return p.size() == n; })) {
            throw std::invalid_argument("The input data (points) is inconsistent: all points must have the same "
                                        "dimension, while I detect differences.");
        }
        // 4 - We check that all l in labels have the same size
        if (!std::all_of(m_labels.begin(), m_labels.end(),
                         [m](const std::vector<double> &l) { return l.size() == m; })) {
            throw std::invalid_argument("The labels are inconsistent: all labels must have the same "
                                        "dimension, while I detect differences.");
        }
        if (m_c == 0) throw std::invalid_argument("Number of columns is 0");
        if (m_r == 0) throw std::invalid_argument("Number of rows is 0");
        if (m_l == 0) throw std::invalid_argument("Number of level-backs is 0");
        if (m_arity < 2) throw std::invalid_argument("Arity must me at least 2.");
        if (m_f.size() == 0) throw std::invalid_argument("Number of basis functions is 0");
    }

    std::vector<std::vector<double>> m_points;
    std::vector<std::vector<double>> m_labels;
    std::vector<std::vector<gdual_d>> m_dpoints;
    std::vector<std::vector<gdual_d>> m_dlabels;
    std::vector<std::string> m_deph_symb;
    unsigned m_r;
    unsigned m_c;
    unsigned m_l;
    unsigned m_arity;
    std::vector<kernel<double>> m_f;
    unsigned m_n_eph;
    unsigned m_parallel_batches;
    // The fact that this is mutable may hamper the performances of the bfe as the thread safetly level
    // of the UDP in pagmo will force copies of the UDP to be made in all threads. This can in principle be
    // avoided, but likely resulting in prepature optimization. (see https://github.com/darioizzo/dcgp/pull/42)
    mutable expression<double> m_cgp;
    mutable expression<gdual_d> m_dcgp;
    mutable std::pair<pagmo::vector_double, pagmo::vector_double> m_cache;
};
} // namespace dcgp
#endif
