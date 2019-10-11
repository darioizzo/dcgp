#ifndef DCGP_GD4CGP_H
#define DCGP_GD4CGP_H

#include <audi/gdual.hpp>
#include <numeric>
#include <pagmo/algorithms/not_population_based.hpp>
#include <pagmo/detail/custom_comparisons.hpp>

#include <pagmo/io.hpp>
#include <pagmo/population.hpp>
#include <stdexcept>

#include <vector>

#include <dcgp/kernel.hpp>
#include <dcgp/kernel_set.hpp>
#include <dcgp/problems/symbolic_regression.hpp>

namespace dcgp
{
class gd4cgp : public pagmo::not_population_based
{
public:
    /// Single entry of the log (iter, fevals, gevals, gradient_magnitude, lr, best)
    typedef std::tuple<unsigned, unsigned long long, unsigned long long, double, double, double> log_line_type;
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
    gd4cgp(unsigned max_iter = 1u, double lr = 1., double lr_min = 1e-3)
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
        auto n_obj = prob.get_nobj();
        const auto bounds = prob.get_bounds();
        auto NP = pop.size();
        auto fevals0 = prob.get_fevals(); // fevals already made
        auto gevals0 = prob.get_gevals(); // gevals already made
        auto count = 1u;                  // regulates the screen output

        // PREAMBLE-------------------------------------------------------------------------------------------------
        // Check whether the problem is suitable for gd4cgp
        // If the UDP in pop is not a symbolic_regression UDP, udp_ptr will be NULL
        auto udp_ptr = prob.extract<symbolic_regression>();
        if (!udp_ptr) {
            throw std::invalid_argument(prob.get_name() + " does not seem to be a symbolic regression problem. "
                                        + get_name()
                                        + " can only be used on problems of the type dcgp::symbolic_regression.");
        }
        if (udp_ptr->get_cgp().get_eph_val().size() == 0u) {
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
            if (m_verbosity > 0u) {
                pagmo::print("Exit condition -- maximum number of iteration is zero.", '\n');
            }
            return pop;
        }
        // ---------------------------------------------------------------------------------------------------------
        // No throws, all valid: we clear the logs
        m_log.clear();

        // 1 - We select the chromosome in the population.
        auto sel_xf = select_individual(pop);
        pagmo::vector_double x0(std::move(sel_xf.first)), fit0(std::move(sel_xf.second));
        pagmo::vector_double x1(x0);
        pagmo::vector_double fit1(fit0);

        // Get out if there is nothing to do.
        if (!std::isfinite(fit0[0])) {
            if (m_verbosity > 0u) {
                pagmo::print("Exit condition -- population best is not finite: ", fit0[0], '\n');
            }
            return pop;
        }

        // 2 - Gradient Descent iterations
        double lr = m_lr;
        double loss_gradient_norm = 0.;

        for (unsigned iter = 1u; iter <= m_max_iter; ++iter) {
            // We log
            if (m_verbosity > 0u) {
                // Every m_verbosity generations print a log line
                if (iter % m_verbosity == 1u || m_verbosity == 1u) {
                    // Every 50 lines print the column names
                    if (count % 50u == 1u) {
                        pagmo::print("\n", std::setw(7), "Iter:", std::setw(15), "Fevals:", std::setw(15),
                                     "Gevals:", std::setw(15), "grad norm:", std::setw(15), "lr:", std::setw(15),
                                     "Best:\n");
                    }
                    log_single_line(iter-1, prob, fevals0, gevals0, loss_gradient_norm, lr, fit0);
                    ++count;
                }
            }
            auto grad = prob.gradient(x0);
            loss_gradient_norm = std::sqrt(std::inner_product(grad.begin(), grad.end(), grad.begin(), 0.));
            if (loss_gradient_norm < 1e-13 || !std::isfinite(loss_gradient_norm)) { // nothing to do for GD
                if (m_verbosity > 0u) {
                    log_single_line(iter, prob, fevals0, gevals0, loss_gradient_norm, lr, fit0);
                    pagmo::print("Exit condition -- vanishing or nan/inf gradient = ", loss_gradient_norm, '\n');
                }
                replace_individual(pop, x0, fit0);
                return pop;
            }
            // The SGD update rule
            std::transform(grad.begin(), grad.end(), x1.data(), x1.data(),
                           [lr, loss_gradient_norm](double a, double b) { return b - a / loss_gradient_norm * lr; });
            fit1 = prob.fitness(x1);
            if (pagmo::detail::less_than_f(fit1[0], fit0[0])) {
                x0 = x1;
                fit0 = fit1;
                lr = lr * 1.5;
            } else {
                lr = lr / 4.;
                x1 = x0;
                if (lr < m_lr_min) {
                    if (m_verbosity > 0u) {
                        log_single_line(iter, prob, fevals0, gevals0, loss_gradient_norm, lr, fit0);
                        pagmo::print("Exit condition -- vanishing learning rate = ", lr, '\n');
                    }
                    replace_individual(pop, x0, fit0);
                    return pop;
                }
            }
        }
        if (m_verbosity > 0u) {
            log_single_line(m_max_iter, prob, fevals0, gevals0, loss_gradient_norm, lr, fit0);
            pagmo::print("Exit condition -- max iterations = ", m_max_iter, '\n');
        }
        replace_individual(pop, x0, fit0);
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
    // This prints to screen and logs one single line.
    inline void log_single_line(unsigned iter, const pagmo::problem &prob, unsigned long long fevals0,
                                unsigned long long gevals0, double loss_gradient_norm, double lr,
                                pagmo::vector_double fit0) const
    {
        pagmo::print(std::setw(7), iter, std::setw(15), prob.get_fevals() - fevals0, std::setw(15),
                     prob.get_gevals() - gevals0, std::setw(15), loss_gradient_norm, std::setw(15), lr, std::setw(15),
                     fit0[0], '\n');
        m_log.emplace_back(iter - 1, prob.get_fevals() - fevals0, prob.get_gevals() - gevals0, loss_gradient_norm, lr,
                           fit0[0]);
    }

    unsigned m_max_iter;
    double m_lr;
    double m_lr_min;
    unsigned m_verbosity;
    mutable log_type m_log;
};
} // namespace dcgp
#endif