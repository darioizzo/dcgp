#ifndef DCGP_GD4CGP_H
#define DCGP_GD4CGP_H

#include <audi/gdual.hpp>
#include <numeric>
#include <pagmo/io.hpp>
#include <pagmo/population.hpp>
#include <stdexcept>
#include <vector>

#include <dcgp/kernel.hpp>
#include <dcgp/kernel_set.hpp>
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
            return pop;
        }
        // ---------------------------------------------------------------------------------------------------------
        // No throws, all valid: we clear the logs
        m_log.clear();

        // 1 - We make a copy of the cgp which we will use to create the gdual version
        auto cgp = udp_ptr->get_cgp();
        unsigned n_eph = static_cast<unsigned>(cgp.get_eph_val().size());
        auto n_in = cgp.get_n() - n_eph;
        auto n_out = cgp.get_m();
        auto n_rows = cgp.get_r();
        auto n_cols = cgp.get_c();
        auto n_lb = cgp.get_l();
        auto arity = cgp.get_arity();
        auto eph_symb = cgp.get_eph_symb();
        std::vector<kernel<double>> kernels_double = cgp.get_f();
        kernel_set<audi::gdual_d> kernel_set;
        for (const auto &f : kernels_double) {
            kernel_set.push_back(f.get_name());
        }
        // The gdual version of the cgp (this will compute over the gdual algebra)
        expression<gdual_d> dcgp(n_in, n_out, n_rows, n_cols, n_lb, arity, kernel_set(), n_eph);
        dcgp.set_eph_symb(eph_symb); // for consistency

        // 2 - We extract the best chromosome in the population.
        auto best_idx = pop.best_idx();
        auto best_x = pop.get_x()[best_idx];
        double best_f = pop.get_f()[best_idx][0];
        // ... and set the dCGP with it
        std::vector<unsigned> best_xu(best_x.size() - n_eph);
        std::transform(best_x.data() + n_eph, best_x.data() + best_x.size(), best_xu.begin(),
                       [](double a) { return boost::numeric_cast<unsigned>(a); });
        dcgp.set(best_xu);
        std::vector<audi::gdual_d> eph_val;
        for (decltype(n_eph) i = 0u; i < n_eph; ++i) {
            eph_val.emplace_back(best_x[i], eph_symb[i], 1u); // Only first derivative is needed
        }
        dcgp.set_eph_val(eph_val);

        // 3 - We Transform the points and labels into gduals as to be able to call loss
        auto points = udp_ptr->get_points();
        auto labels = udp_ptr->get_labels();
        auto parallel_batches = udp_ptr->get_parallel_batches();
        std::vector<std::vector<audi::gdual_d>> points_gdual;
        std::vector<std::vector<audi::gdual_d>> labels_gdual;
        for (const auto &point : points) {
            std::vector<audi::gdual_d> point_gdual;
            for (const auto &item : point) {
                point_gdual.emplace_back(item);
            }
            points_gdual.push_back(point_gdual);
        }
        for (const auto &label : labels) {
            std::vector<audi::gdual_d> label_gdual;
            for (const auto &item : label) {
                label_gdual.emplace_back(item);
            }
            labels_gdual.push_back(label_gdual);
        }

        // We create the symbol set of the differentials
        std::vector<std::string> eph_symb_d; // symbol set needs to be with d added
        for (const auto &symb : eph_symb) {
            eph_symb_d.push_back("d" + symb);
        }

        // Gradient Descent iterations
        double lr = m_lr;
        gdual_d loss;

        for (unsigned i = 0u; i < 100; ++i) {
            loss = dcgp.loss(points_gdual, labels_gdual, "MSE", parallel_batches);
            loss.extend_symbol_set(eph_symb_d); // when the eph variables are not in the expression we still need their
                                                // symbols as to get a derivative
            audi::print(i, " ", loss.constant_cf(), "\n");
            std::vector<double> loss_gradient(n_eph, 0.);
            if (!(loss.get_order() == 0)) { // this happens when input terminals of the eph constants are inactive
                                            // (gradient is then zero)
                for (decltype(n_eph) i = 0u; i < n_eph; ++i) {
                    std::vector<double> coeff(n_eph, 0.);
                    coeff[i] = 1.;
                    loss_gradient[i] = loss.get_derivative(coeff);
                }
            }
            double loss_gradient_norm
                = std::sqrt(std::inner_product(loss_gradient.begin(), loss_gradient.end(), loss_gradient.begin(), 0.));
            // audi::print(loss_gradient, loss_gradient_norm);
            if (loss_gradient_norm == 0. || !std::isfinite(loss_gradient_norm)) { // nothing to do for GD
                return pop;
            }
            // The SGD update rule
            std::transform(loss_gradient.begin(), loss_gradient.end(), eph_val.begin(), eph_val.begin(),
                           [lr, loss_gradient_norm](double a, gdual_d &b) { return b - a / loss_gradient_norm * lr; });

            dcgp.set_eph_val(eph_val);
        }

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