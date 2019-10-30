#ifndef DCGP_SYMBOLIC_REGRESSION_H
#define DCGP_SYMBOLIC_REGRESSION_H
#include <audi/gdual.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/range/algorithm/transform.hpp>
#include <numeric> // std::accumulate
#include <pagmo/io.hpp>
#include <pagmo/population.hpp>
#include <pagmo/types.hpp>
#include <symengine/expression.h>
#include <vector>

#include <dcgp/expression.hpp>
#include <dcgp/kernel_set.hpp>
#include <dcgp/rng.hpp>

namespace dcgp
{
/// A Symbolic Regression problem
/**
 *
 * Symbolic regression is a type of regression analysis that searches the space of mathematical expressions to 
 * find the model that best fits a given dataset, both in terms of accuracy and simplicity 
 * (ref: https://en.wikipedia.org/wiki/Symbolic_regression). It also is one of the core applications
 * for Differentiable Cartesian Genetic Programming.
 *
 * This class provides an easy way to instantiate symbolic regression problems as optimization problems having
 * a continuous part (i.e. the value of the parameters in the model) and an integer part (i.e. the representation of
 * the model computational graph). The instantiated object can be used as UDP (User Defined Problem) in the pagmo optimization suite.
 *
 * The symbolic regression problem can be instantiated both as a single and a two-objectives problem. In the second
 * case, aside the Mean Squared Error, the formula complexity will be considered as an objective.
 *
 */
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
        : m_points(1), m_labels(1), m_r(1), m_c(1), m_l(1), m_arity(2), m_f(kernel_set<double>({"sum"})()), m_n_eph(0),
          m_multi_objective(true), m_parallel_batches(0u)
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
                        unsigned r = 1u,     // n. rows
                        unsigned c = 10u,    // n. columns
                        unsigned l = 11u,    // n. levels-back
                        unsigned arity = 2u, // basis functions' arity
                        std::vector<kernel<double>> f
                        = kernel_set<double>({"sum", "diff", "mul", "pdiv"})(), // functions
                        unsigned n_eph = 0u,                                    // number of ephemeral constants
                        bool multi_objective = false,  // when true the fitness also returns the formula complexity
                        unsigned parallel_batches = 0u // number of parallel batches
                        )
        : m_points(points), m_labels(labels), m_r(r), m_c(c), m_l(l), m_arity(arity), m_f(f), m_n_eph(n_eph),
          m_multi_objective(multi_objective), m_parallel_batches(parallel_batches)
    {
        unsigned n;
        unsigned m;
        // We check the input against obvious sanity criteria.
        sanity_checks(n, m);
        // We initialize the inner cgp expression
        auto seed = random_device::next();
        m_cgp = expression<double>(n, m, m_r, m_c, m_l, m_arity, m_f, m_n_eph, seed);
        // We initialize the inner dcgp expression
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
        // We create the symbols of the input variables here
        for (decltype(m_points[0].size()) i = 0u; i < m_points[0].size(); ++i) {
            m_symbols.push_back("x" + std::to_string(i));
        }
    }

    /// Number of objectives
    /**
     * Returns the number of objectives.
     *
     * @return the number of objectives.
     */
    pagmo::vector_double::size_type get_nobj() const
    {
        return 1u + m_multi_objective;
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
        std::vector<double> retval(1u + m_multi_objective, 0);
        // Here we set the CGP member from the chromosome
        set_cgp(x);
        if (x == m_cache_fitness.first) {
            retval[0] = m_cache_fitness.second[0];
        } else {
            // And we compute the MSE loss splitting the data in n batches.
            // TODO: make this work also when m_parallel_batches does not divide exactly the data size.
            retval[0] = m_cgp.loss(m_points, m_labels, "MSE", m_parallel_batches);
        }
        // In the multiobjective case we compute the formula complexity
        if (m_multi_objective) {
            std::ostringstream ss;
            // A first "naive" implementation of the formula complexity measure is the length of
            // the shortest string among pretty and prettier. That is among the raw cgp expression
            // and the result of constructing a symengine expression out of it (which carries out some
            // basic simplifications but that may results in rare occasions in a longer string).
            std::vector<std::string> pretty = m_cgp(m_symbols);
            double l_pretty = std::accumulate(pretty.begin(), pretty.end(), 0., [](double a, std::string b) {
                return a + static_cast<double>(b.length());
            });
            double l_prettier = 0.;
            for (decltype(pretty.size()) i = 0u; i < pretty.size(); ++i) {
                SymEngine::Expression prettier(pretty[i]);
                pagmo::stream(ss, prettier);
                auto string = ss.str();
                // We remove whitespaces too
                l_prettier += static_cast<double>(string.length() - static_cast<decltype(string.length())>(std::count(string.begin(), string.end(), ' ')));
            }
            // Here we define the formula complexity
            retval[1] = std::min(l_pretty, l_prettier);
        }
        return retval;
    }

    /// Gradient computation
    /**
     * Computes the gradient of the loss with respect to the ephemeral constants (i.e. the continuous part of the
     * chromosome).
     *
     * @param x the decision vector.
     *
     * @return the gradient in \p x.
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
        // Now we extract fitness and gradient and store the values in th caches or in the return value
        m_cache_fitness.first = x;
        m_cache_fitness.second.resize(get_nobj());
        m_cache_fitness.second[0] = loss.constant_cf();
        loss.extend_symbol_set(m_deph_symb);
        if (!(loss.get_order() == 0u)) { // this happens when input terminals of the eph constants are inactive
                                         // (gradient is then zero)
            for (decltype(m_n_eph) i = 0u; i < m_n_eph; ++i) {
                std::vector<unsigned> coeff(m_n_eph, 0.);
                coeff[i] = 1.;
                retval[i] = loss.get_derivative(coeff);
            }
        }
        return retval;
    }

    /// Sparsity pattern (gradient)
    /**
     * Returns the sparsity pattern of the gradient. The sparsity patter is dense in the continuous part of the
     * chromosome. (this is a result of assuming all ephemeral constants are actually in the expression, if not zeros
     * will be returned)
     *
     * @return the gradient sparsity pattern.
     */
    pagmo::sparsity_pattern gradient_sparsity() const
    {
        pagmo::sparsity_pattern retval;
        for (decltype(m_n_eph) i = 0u; i < m_n_eph; ++i) {
            retval.push_back({0u, i});
        }
        return retval;
    }

    /// Hessian computation
    /**
     * Computes the hessian of the loss with respect to the ephemeral constants (i.e. the continuous part of the
     * chromosome).
     *
     * @param x the decision vector.
     *
     * @return the hessian in \p x.
     */
    std::vector<pagmo::vector_double> hessians(const pagmo::vector_double &x) const
    {
        // The hessian sparsity and dimension
        auto hs = hessians_sparsity();
        auto hd = hs[0].size();
        // Initializing the return value (hessian) to zeros.
        std::vector<pagmo::vector_double> retval;
        for (const auto &item : hs) {
            retval.emplace_back(item.size(), 0.);
        }
        // Initializing the gradient to zeros.
        m_cache_gradient.first = x;
        m_cache_gradient.second.clear();
        m_cache_gradient.second.resize(m_n_eph, 0.);

        // The chromosome has a floating point part (the ephemeral constants) and an integer part (the encoded CGP).
        // 1 - We extract the integer part and represent it as an unsigned vector to set the CGP expression.
        std::vector<unsigned> xu(x.size() - m_n_eph);
        std::transform(x.data() + m_n_eph, x.data() + x.size(), xu.begin(),
                       [](double a) { return boost::numeric_cast<unsigned>(a); });
        m_dcgp.set(xu);
        // 2 - We use the floating point part of the chromosome to set the ephemeral constants values of the dcgp
        // member.
        std::vector<audi::gdual_d> eph_val;
        for (decltype(m_n_eph) i = 0u; i < m_n_eph; ++i) {
            eph_val.emplace_back(x[i], m_dcgp.get_eph_symb()[i], 2u); // First and second order derivative are needed
        }
        m_dcgp.set_eph_val(eph_val);
        // 3 - We compute the MSE loss and its differentials.
        auto loss = m_dcgp.loss(m_dpoints, m_dlabels, "MSE", m_parallel_batches);
        // We make sure all symbols are in so that we get zeros when querying for a variable not in the gdual
        loss.extend_symbol_set(m_deph_symb);

        // Now we extract fitness, gradient and hessians from the gdual and store the values (retval or cache)
        m_cache_fitness.first = x;
        m_cache_fitness.second.resize(get_nobj());
        m_cache_fitness.second[0] = loss.constant_cf();
        // We compute the gradient and the hessian only if
        // the loss depends on at least one ephemeral constant.
        // Otherwise the initialization values will be returned, that is zeros.
        if (!(loss.get_order() == 0u)) {
            // gradient (we cache it)
            for (decltype(m_n_eph) i = 0u; i < m_n_eph; ++i) {
                std::vector<unsigned> coeff(m_n_eph, 0.);
                coeff[i] = 1.;
                m_cache_gradient.second[i] = loss.get_derivative(coeff);
            }
            // hessian (we return it)
            for (decltype(hd) i = 0u; i < hd; ++i) {
                std::vector<unsigned> coeff(m_n_eph, 0.);
                coeff[hs[0][i].first] = 1;
                coeff[hs[0][i].second] += 1;
                retval[0][i] = loss.get_derivative(coeff);
            }
        }
        return retval;
    }

    /// Sparsity pattern (hessian)
    /**
     * Returns the sparsity pattern of the hessian. The sparsity patter is dense in the continuous part of the
     * chromosome. (this is a result of assuming all ephemeral constants are actually in the expression, if not zeros
     * will be returned)
     *
     * @return the hessian sparsity pattern.
     */
    std::vector<pagmo::sparsity_pattern> hessians_sparsity() const
    {
        std::vector<pagmo::sparsity_pattern> retval;
        pagmo::sparsity_pattern sp;
        for (decltype(m_n_eph) i = 0u; i < m_n_eph; ++i) {
            for (decltype(i) j = 0u; j <= i; ++j) {
                sp.push_back({i, j});
            }
        }
        retval.push_back(std::move(sp));
        if (m_multi_objective) {
            retval.push_back({{}});
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
     * Returns the integer dimension of the problem.
     *
     * @return the integer dimension of the problem.
     */
    pagmo::vector_double::size_type get_nix() const
    {
        return m_cgp.get_lb().size();
    }

    /// Problem name
    /**
     * Returns a string containing the problem name.
     *
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
        pagmo::stream(ss, "\tData dimension (in): ", m_points[0].size(), "\n");
        pagmo::stream(ss, "\tData dimension (out): ", m_labels[0].size(), "\n");
        pagmo::stream(ss, "\tData size: ", m_points.size(), "\n");
        pagmo::stream(ss, "\tKernels: ", m_cgp.get_f(), "\n");
        return ss.str();
    }

    /// Human-readable representation of a decision vector.
    /**
     * A human readable representation of the chromosome is here obtained by calling directly
     * the expression::operator() assuming as inputs variables names \f$x_1, x_2, ...\f$ and
     * as ephemeral constants names \f$c_1, c_2, ...\f$
     *
     * @param[in] x a valid chromosome.
     *
     * @return a string containing the mathematical expression represented by *x*.
     */
    std::string pretty(const pagmo::vector_double &x) const
    {
        // Here we set the CGP data member from the chromosome
        set_cgp(x);

        std::ostringstream ss;
        std::vector<std::string> symbols;
        for (decltype(m_points[0].size()) i = 0u; i < m_points[0].size(); ++i) {
            symbols.push_back("x" + std::to_string(i));
        }
        pagmo::stream(ss, m_cgp(symbols));
        return ss.str();
    }

    /// Human-readable representation of a decision vector.
    /**
     * A human readable representation of the chromosome is here obtained by using symengine
     * the expression::operator() assuming as inputs variables names \f$x_1, x_2, ...\f$ and
     * as ephemeral constants names \f$c_1, c_2, ...\f$
     *
     * @param[in] x a valid chromosome.
     *
     * @return a string containing the mathematical expression represented by *x*.
     */
    std::string prettier(const pagmo::vector_double &x) const
    {
        // Here we set the CGP data member from the chromosome
        set_cgp(x);

        std::ostringstream ss;
        std::vector<std::string> symbols;
        for (decltype(m_points[0].size()) i = 0u; i < m_points[0].size(); ++i) {
            symbols.push_back("x" + std::to_string(i));
        }
        auto raws = m_cgp(symbols);
        std::vector<SymEngine::Expression> exs;
        for (auto const &raw : raws) {
            exs.emplace_back(raw);
        }
        pagmo::stream(ss, exs);
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

    pagmo::thread_safety get_thread_safety() const
    {
        return pagmo::thread_safety::none;
    }

private:
    // This setter can be marked const as m_cgp is mutable
    void set_cgp(const pagmo::vector_double &x) const
    {
        // We need to make a copy of the chromosome as to represents its genes as unsigned
        std::vector<unsigned> xu(x.size() - m_n_eph);
        std::transform(x.data() + m_n_eph, x.data() + x.size(), xu.data(),
                       [](double a) { return boost::numeric_cast<unsigned>(a); });
        m_cgp.set(xu);
        // 2 - We set the floating point part as ephemeral constants.
        std::vector<double> eph_val(x.data(), x.data() + m_n_eph);
        m_cgp.set_eph_val(eph_val);
    }

    void sanity_checks(unsigned &n, unsigned &m) const
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
    std::vector<std::string> m_symbols;

    unsigned m_r;
    unsigned m_c;
    unsigned m_l;
    unsigned m_arity;
    std::vector<kernel<double>> m_f;
    unsigned m_n_eph;
    bool m_multi_objective;
    unsigned m_parallel_batches;
    // The fact that this is mutable may hamper the performances of the bfe as the thread safetly level
    // of the UDP in pagmo will force copies of the UDP to be made in all threads. This can in principle be
    // avoided, but likely resulting in prepature optimization. (see https://github.com/darioizzo/dcgp/pull/42)
    mutable expression<double> m_cgp;
    // TODO: this should be vectorized gduals
    mutable expression<gdual_d> m_dcgp;
    mutable std::pair<pagmo::vector_double, pagmo::vector_double> m_cache_fitness;
    mutable std::pair<pagmo::vector_double, pagmo::vector_double> m_cache_gradient;
}; // namespace dcgp
} // namespace dcgp
#endif
