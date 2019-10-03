#ifndef DCGP_GD4CGP_H
#define DCGP_GD4CGP_H

#include <stdexcept>
#include <vector>
#include <pagmo/population.hpp>
#include <pagmo/io.hpp>

#include <dcgp/problems/symbolic_regression.hpp>

namespace dcgp
{
class gd4cgp
{
public:
    /// Single entry of the log (iter, fevals, best, step_size)
    typedef std::tuple<unsigned, unsigned long long, double, double> log_line_type;
    /// The log
    typedef std::vector<log_line_type> log_type;

    /// Constructor
    /**
     * Constructs a stochastic gradient descent for use with a cgp::symbolic_regression UDP that has
     * any strictly positive number of ephemeral constants
     *
     * @param max_iter maximum number of iterations.
     * @param lr initial learning rate.
     * @param lr_min minimum learning rate.
     *
     * @throws std::invalid_argument if *lr_min* is not in (0, *lr*)
     */
    gd4cgp(unsigned max_iter = 1u, double lr = 0.1, double lr_min = 1e-3)
        : m_max_iter(max_iter), m_lr(lr), m_lr_min(lr_min), m_verbosity(0u)
    {
        if (lr_min <= 0.) {
            throw std::invalid_argument("The minimum learning rate must be strictly positive.");
        }
        if (lr_min >= lr) {
            throw std::invalid_argument(
                "The minimum learning rate must be strictly smaller than the initial learning rate.");
        }
    }
    // Algorithm evolve method
    pagmo::population evolve(pagmo::population pop) const
    {
        const auto &prob = pop.get_problem();
        auto dim = prob.get_nx();
        auto n_obj = prob.get_nobj();
        const auto bounds = prob.get_bounds();
        auto NP = pop.size();
        auto fevals0 = prob.get_fevals(); // fevals already made
        auto count = 1u;                  // regulates the screen output

        // PREAMBLE-------------------------------------------------------------------------------------------------
        // Check whether the problem is suitable for es4cgp
        // If the UDP in pop is not a symbolic_regression UDP, udp_ptr will be NULL
        auto udp_ptr = prob.extract<symbolic_regression>();
        if (!udp_ptr) {
            throw std::invalid_argument(prob.get_name() + " does not seem to be a symbolic regression problem. "
                                        + get_name()
                                        + " can only be used on problems of the type dcgp::symbolic_regression.");
        }
        if (udp_ptr->get_eph_val().size() == 0u) {
            throw std::invalid_argument(prob.get_name() + " does not seem to have any ephemeral constants. "
                                        + get_name()
                                        + " only acts on ephemeral constants and needs at least one to be present.");
        }
        if (n_obj > 1) {
            throw std::invalid_argument(prob.get_name() + " has multiple objectives. " + get_name()
                                        + " can only be used on problems that are single objective.");
        }
        if (NP < 1u) {
            throw std::invalid_argument(get_name() + " needs at least 2 individuals in the population, "
                                        + std::to_string(NP) + " detected.");
        }
        // Get out if there is nothing to do.
        if (m_max_iter == 0u) {
            return pop;
        }
        // ---------------------------------------------------------------------------------------------------------

        // No throws, all valid: we clear the logs
        m_log.clear();
        // We make a copy of the cgp which we will use to make mutations.
        auto cgp = udp_ptr->get_cgp();
        // We get the best chromosome in the population.
        auto best_idx = pop.best_idx();
        auto best_x = pop.get_x()[best_idx];
        double best_f = pop.get_f()[best_idx][0];
        // ... and we transform it into unsigned
        std::vector<unsigned> best_xu(best_x.size());
        std::transform(best_x.begin(), best_x.end(), best_xu.begin(),
                       [](double a) { return boost::numeric_cast<unsigned>(a); });
        return pop;
    }

    /// Sets the algorithm verbosity
    /**
     * Sets the verbosity level of the screen output and of the
     * log returned by get_log(). \p level can be:
     * - 0: no verbosity
     * - >0: will print and log one line each \p level generations.
     *
     * Example (verbosity 100):
     * @code{.unparsed}
     *  Gen:        Fevals:         Best:
     *     1              0       0.922363
     *     2              4       0.153203
     *     3              8       0.125378
     *     4             12       0.125378
     *     5             16       0.125378
     *     6             20       0.125378
     * @endcode
     * Gen is the generation number, Fevals the number of function evaluation used, Best is the best fitness found.
     *
     * @param level verbosity level
     */
    void set_verbosity(unsigned level)
    {
        m_verbosity = level;
    }

    /// Gets the verbosity level
    /**
     * @return the verbosity level
     */
    unsigned get_verbosity() const
    {
        return m_verbosity;
    }

    /// Algorithm name
    /**
     * @return a string containing the algorithm name
     */
    std::string get_name() const
    {
        return "GD for CGP: gradient descent for Cartesian Genetic Programming";
    }

    // Extra info
    std::string get_extra_info() const
    {
        std::ostringstream ss;
        pagmo::stream(ss, "\tMaximum number of iterations: ", m_max_iter);
        pagmo::stream(ss, "\n\tInitial learning rate: ", m_lr);
        pagmo::stream(ss, "\n\tMinimum learning rate: ", m_lr_min);
        pagmo::stream(ss, "\n\tVerbosity: ", m_verbosity);
        return ss.str();
    }

    /// Get log
    /**
     * A log containing relevant quantities monitoring the last call to evolve. Each element of the returned
     * <tt>std::vector</tt> is a bee_colony::log_line_type containing: Gen, Fevals, Current best, Best as
     * described in bee_colony::set_verbosity().
     *
     * @return an <tt> std::vector</tt> of bee_colony::log_line_type containing the logged values Gen, Fevals, Current
     * best, Best
     */
    const log_type &get_log() const
    {
        return m_log;
    }

private:
    unsigned m_max_iter;
    double m_lr;
    double m_lr_min;
    unsigned m_verbosity;
    mutable log_type m_log;
};
} // namespace dcgp
#endif