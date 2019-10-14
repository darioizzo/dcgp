#ifndef DCGP_ES4CGP_H
#define DCGP_ES4CGP_H

#include <pagmo/bfe.hpp>
#include <pagmo/detail/custom_comparisons.hpp>
#include <pagmo/io.hpp>
#include <pagmo/population.hpp>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include <dcgp/problems/symbolic_regression.hpp>
#include <dcgp/rng.hpp>

namespace dcgp
{
class es4cgp
{
public:
    /// Single entry of the log (gen, fevals, best, formula)
    typedef std::tuple<unsigned, unsigned long long, double, std::string> log_line_type;
    /// The log
    typedef std::vector<log_line_type> log_type;

    /// Constructor
    /**
     * Constructs an evolutionary strategy algorithm for use with a cgp::symbolic_regression UDP.
     *
     * @param gen number of generations.
     * @param mut_n number of active genes to be mutated.
     * @param ftol the algorithm will exit when the loss is below this tolerance.
     * @param seed seed used by the internal random number generator (default is random).
     *
     * @throws std::invalid_argument if *mut_n* is 0 or *ftol* is negative
     */
    es4cgp(unsigned gen = 1u, unsigned mut_n = 2u, double ftol = 1e-4, unsigned seed = random_device::next())
        : m_gen(gen), m_mut_n(mut_n), m_ftol(ftol), m_e(seed), m_seed(seed), m_verbosity(0u)
    {
        if (mut_n == 0u) {
            throw std::invalid_argument("The number of active mutations is zero, it must be at least 1.");
        }
        if (ftol < 0.) {
            throw std::invalid_argument("The ftol is negative, it must be positive or zero.");
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
        auto udp_ptr = prob.extract<symbolic_regression>();
        // PREAMBLE-------------------------------------------------------------------------------------------------
        // Check whether the problem is suitable for es4cgp
        // If the UDP in pop is not a symbolic_regression UDP, udp_ptr will be NULL
        if (!udp_ptr) {
            throw std::invalid_argument(prob.get_name() + " does not seem to be a symbolic regression problem. "
                                        + get_name()
                                        + " can only be used on problems of the type dcgp::symbolic_regression ");
        }
        if (n_obj > 1) {
            throw std::invalid_argument(prob.get_name() + " has multiple objectives. " + get_name()
                                        + " can only be used on problems that are single objective.");
        }
        if (NP < 2u) {
            throw std::invalid_argument(get_name() + " needs at least 2 individuals in the population, "
                                        + std::to_string(NP) + " detected");
        }
        // Get out if there is nothing to do.
        if (m_gen == 0u) {
            return pop;
        }
        // ---------------------------------------------------------------------------------------------------------

        // No throws, all valid: we clear the logs
        m_log.clear();
        // We make a copy of the cgp which we will use to make mutations.
        auto cgp = udp_ptr->get_cgp();
        // How many ephemeral constants? And their values?
        auto eph_val = cgp.get_eph_val();
        auto n_eph = eph_val.size();
        // We get the best chromosome in the population.
        auto best_idx = pop.best_idx();
        auto best_x = pop.get_x()[best_idx];
        double best_f = pop.get_f()[best_idx][0];
        // ... and we transform it into unsigned
        std::vector<unsigned> best_xu(best_x.size() - n_eph);
        std::transform(best_x.data() + n_eph, best_x.data() + best_x.size(), best_xu.begin(),
                       [](double a) { return boost::numeric_cast<unsigned>(a); });
        // A contiguous vector of chromosomes/fitness vectors for pagmo::bfe input\output is allocated here.
        pagmo::vector_double dvs(NP * dim);
        pagmo::vector_double fs(NP * n_obj);
        // We fill in the chromosomes the ephemeral constant values (they will not change anymore during this evolve)
        for (decltype(NP) i = 0u; i < NP; ++i) {
            std::copy(eph_val.begin(), eph_val.end(), dvs.data() + dim * i);
        }

        // Main loop
        for (decltype(m_gen) gen = 1u; gen <= m_gen; ++gen) {
            // Logs and prints (verbosity modes > 1: a line is added every m_verbosity generations)
            if (m_verbosity > 0u) {
                // Every m_verbosity generations print a log line
                if (gen % m_verbosity == 1u || m_verbosity == 1u) {
                    // Every 50 lines print the column names
                    if (count % 50u == 1u) {
                        pagmo::print("\n", std::setw(7), "Gen:", std::setw(15), "Fevals:", std::setw(15),
                                     "Best:", "\tFormula:\n");
                    }
                    auto formula = udp_ptr->prettier(best_x);
                    log_single_line(gen - 1, prob.get_fevals() - fevals0, best_f, formula);
                    ++count;
                }
            }

            // 1 - We generate new NP individuals mutating the best and we write on the dvs for pagmo::bfe to evaluate
            // their fitnesses.
            for (decltype(NP) i = 0u; i < NP; ++i) {
                // To mutate active chromosomes we use a copy of the cgp of the UDP.
                cgp.set(best_xu);
                cgp.mutate_active(m_mut_n); // TODO: check if this is good
                // cgp.mutate(m_mut_n + i);
                std::vector<unsigned> mutated_x = cgp.get();
                std::transform(mutated_x.begin(), mutated_x.end(), dvs.data() + i * dim + n_eph,
                               [](unsigned a) { return boost::numeric_cast<double>(a); });
            }
            // 2 - We compute their fitnesses calling the bfe
            fs = m_bfe(prob, dvs);
            // Check if ftol exit condition is met
            if (best_f < m_ftol) {
                if (m_verbosity > 0u) {
                    auto formula = udp_ptr->prettier(best_x);
                    log_single_line(gen, prob.get_fevals() - fevals0, best_f, formula);
                    ++count;
                    pagmo::print("Exit condition -- ftol < ", m_ftol, "\n");
                }
                update_pop(pop, dvs, fs, best_f, best_xu, NP, dim, eph_val);
                return pop;
            }
            // 4 - We reinsert the mutated individuals in the population if their fitness is
            // less than, or equal, to the one from the parent.
            for (decltype(NP) i = 0u; i < NP; ++i) {
                if (pagmo::detail::less_than_f(fs[i], best_f)) {
                    best_f = fs[i];
                    std::copy(dvs.data() + i * dim, dvs.data() + (i + 1) * dim, best_x.begin());
                    // best_xu is updated here
                    std::transform(dvs.data() + i * dim + n_eph, dvs.data() + (i + 1) * dim, best_xu.begin(),
                                   [](double a) { return boost::numeric_cast<unsigned>(a); });
                }
            }
        }
        // Evolution has terminated and we now update into the pagmo::pop.
        update_pop(pop, dvs, fs, best_f, best_xu, NP, dim, eph_val);
        // At the end, the pagmo::population will contain the best individual together with its best NP-1 mutants.
        // We log the last iteration
        if (m_verbosity > 0u) {
            auto formula = udp_ptr->prettier(best_x);
            log_single_line(m_gen, prob.get_fevals() - fevals0, best_f, formula);
            pagmo::print("Exit condition -- generations = ", m_gen, '\n');
        }
        return pop;
    }

    // Sets the seed
    void set_seed(unsigned seed)
    {
        m_e.seed(seed);
        m_seed = seed;
    }

    /// Gets the seed
    /**
     * @return the seed controlling the algorithm stochastic behaviour
     */
    unsigned get_seed() const
    {
        return m_seed;
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
        return "ES for CGP: Evolutionary strategy for Cartesian Genetic Programming";
    }

    // Extra info
    std::string get_extra_info() const
    {
        std::ostringstream ss;
        pagmo::stream(ss, "\tMaximum number of generations: ", m_gen);
        pagmo::stream(ss, "\n\tNumber of active mutations: ", m_mut_n);
        pagmo::stream(ss, "\n\tExit condition of the final loss (ftol): ", m_ftol);
        pagmo::stream(ss, "\n\tVerbosity: ", m_verbosity);
        pagmo::stream(ss, "\n\tSeed: ", m_seed);
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
    void log_single_line(unsigned gen, unsigned long long fevals, double best_f, const std::string &formula) const
    {
        pagmo::print(std::setw(7), gen, std::setw(15), fevals, std::setw(15), best_f, "\t",
                     formula.substr(0, 40) + " ...", '\n');
        m_log.emplace_back(gen, fevals, best_f, formula);
    }

    // Used to update the population from the dvs, fs used via the bfe. Typically done at the end of the evolve when an
    // exit condition is met.
    void update_pop(pagmo::population &pop, const pagmo::vector_double &dvs, const pagmo::vector_double &fs,
                    double best_f, const std::vector<unsigned> &best_xu, pagmo::vector_double::size_type NP,
                    pagmo::vector_double::size_type dim, const std::vector<double> &eph_val) const
    {
        pagmo::vector_double best_x(dim);
        // First, the latest generation of mutants (if better)
        for (decltype(NP) i = 0u; i < NP; ++i) {
            if (pagmo::detail::less_than_f(fs[i], pop.get_f()[i][0])) {
                pop.set_xf(i, pagmo::vector_double(dvs.data() + i * dim, dvs.data() + (i + 1) * dim),
                           pagmo::vector_double(1., fs[i]));
            }
        }
        auto worst_idx = pop.worst_idx();
        double worst_f_in_pop = pop.get_f()[worst_idx][0];
        // Then the best in place of the worst
        if (pagmo::detail::less_than_f(best_f, worst_f_in_pop)) {
            std::copy(eph_val.begin(), eph_val.end(), best_x.begin());
            std::transform(best_xu.begin(), best_xu.end(), best_x.data() + eph_val.size(),
                           [](unsigned a) { return boost::numeric_cast<double>(a); });
            pop.set_xf(worst_idx, best_x, pagmo::vector_double(1, best_f));
        }
    }
    unsigned m_gen;
    unsigned m_mut_n;
    double m_ftol;
    mutable detail::random_engine_type m_e;
    unsigned m_seed;
    unsigned m_verbosity;
    mutable log_type m_log;
    pagmo::bfe m_bfe;
};
} // namespace dcgp
#endif
