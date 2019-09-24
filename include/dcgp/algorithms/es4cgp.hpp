#ifndef DCGP_ES4CGP_H
#define DCGP_ES4CGP_H

#include <pagmo/bfe.hpp>
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
    /// Single entry of the log (gen, fevals, best, df)
    typedef std::tuple<unsigned, unsigned long long, double, double> log_line_type;
    /// The log
    typedef std::vector<log_line_type> log_type;

    /// Constructor
    /**
     * Constructs an evolutionary strategy algorithm for use with a cgp::symbolic_regression UDP.
     *
     * @param gen number of generations.
     * @param mutation_type one of "active" or "any". If "active" only active genes will be mutated, otherwise all genes
     * can be mutated.
     * @param mut_n number of active genes to be mutated (when *mutation_type* is "any" this parameter has no effect).
     * @param mut_fraction fraction of the genes to be mutated. The number of genes mutated will be determined rounding
     * up to the nearest integer >= 1 (when *mutation_type* is "active" this parameter has no effect)..
     * @param seed seed used by the internal random number generator (default is random)
     *
     * @throws std::invalid_argument if limit equals 0
     */
    es4cgp(unsigned gen = 1, std::string mutation_type = "active", unsigned mut_n = 2, double mut_fraction = 0.1,
           unsigned seed = random_device::next())
        : m_gen(gen), m_mutation_type(mutation_type), m_mut_n(mut_n), m_mut_fraction(mut_fraction), m_e(seed),
          m_seed(seed), m_verbosity(0u)
    {
        // TODO: add checks
    }
    // Algorithm evolve method
    pagmo::population evolve(pagmo::population pop) const
    {
        const auto &prob = pop.get_problem();
        auto dim = prob.get_nx();
        const auto bounds = prob.get_bounds();
        const auto &lb = bounds.first;
        const auto &ub = bounds.second;
        auto NP = pop.size();
        auto fevals0 = prob.get_fevals(); // fevals already made
        auto count = 1u;                  // regulates the screen output

        // PREAMBLE-------------------------------------------------------------------------------------------------
        // Check whether the problem is suitable for es4cgp
        // If the UDP in pop is not a symbolic_regression UDP, udp_ptr will be NULL
        auto udp_ptr = prob.extract<symbolic_regression>();
        if (!udp_ptr) {
            throw std::invalid_argument(prob.get_name() + " does not seem to be a symbolic regression problem."
                                        + get_name()
                                        + " can only be used on problems of the type dcgp::symbolic_regression ");
        }
        if (NP < 2u) {
            throw std::invalid_argument(prob.get_name() + " needs at least 2 individuals in the population, "
                                        + std::to_string(NP) + " detected");
        }
        // Get out if there is nothing to do.
        if (m_gen == 0u) {
            return pop;
        }
        // ---------------------------------------------------------------------------------------------------------

        // No throws, all valid: we clear the logs
        m_log.clear();
        // Main loop.
        for (decltype(m_gen) gen = 1u; gen <= m_gen; ++gen) {
            // Logs and prints (verbosity modes > 1: a line is added every m_verbosity generations)
            if (m_verbosity > 0u) {
                // Every m_verbosity generations print a log line
                if (gen % m_verbosity == 1u || m_verbosity == 1u) {
                    auto best_idx = pop.best_idx();
                    auto worst_idx = pop.worst_idx();

                    // Every 50 lines print the column names
                    if (count % 50u == 1u) {
                        pagmo::print("\n", std::setw(7), "Gen:", std::setw(15), "Fevals:", std::setw(15),
                                     "Best:", std::setw(15), "Current Best:\n");
                    }
                    pagmo::print(std::setw(7), gen, std::setw(15), prob.get_fevals() - fevals0, std::setw(15),
                                 pop.champion_f()[0], std::setw(15),
                                 pop.get_f()[worst_idx][0] - pop.get_f()[best_idx][0], '\n');
                    ++count;
                    // Logs
                    m_log.emplace_back(gen, prob.get_fevals() - fevals0, pop.champion_f()[0], pop.get_f()[best_idx][0]);
                }
            }
            // 1 - We get the best chromosome.
            auto best_x = pop.get_x()[pop.best_idx()];
            std::vector<unsigned> xu(best_x.size());
            std::transform(best_x.begin(), best_x.end(), xu.begin(),
                           [](double a) { return boost::numeric_cast<unsigned>(a); });
            // 2 - We mutate it and record the contiguous values of the mutated chromosomes for bfe.
            pagmo::vector_double dvs(NP * dim);
            for (decltype(NP) i = 0u; i < NP; ++i) {
                auto cgp = udp_ptr->get_cgp();
                cgp.set(xu);
                cgp.mutate_active(m_mut_n);
                std::vector<unsigned> mutated_x = cgp.get();
                std::transform(mutated_x.begin(), mutated_x.end() , dvs.begin() + i * dim,
                               [](unsigned a) { return boost::numeric_cast<double>(a); });
            }
            // 3 - We compute their fitnesses
            auto fs = m_bfe(prob, dvs);
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
     *     Gen:        Fevals:          Best:  Current Best:
     *        1             40         261363         261363
     *      101           4040        112.237        267.969
     *      201           8040        20.8885        265.122
     *      301          12040        20.6076        20.6076
     *      401          16040         18.252        140.079
     * @endcode
     * Gen is the generation number, Fevals the number of function evaluation used, , Best is the best fitness found,
     * Current best is the best fitness currently in the population.
     *
     * @param level verbosity level
     */
    void set_verbosity(unsigned level)
    {
        m_verbosity = level;
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
        pagmo::stream(ss, "\n\tMutation type: ", m_mutation_type);
        pagmo::stream(ss, "\n\tNumber of active mutations: ", m_mut_n);
        pagmo::stream(ss, "\n\tMutation fraction: ", m_mut_fraction);
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
    unsigned m_gen;
    std::string m_mutation_type;
    unsigned m_mut_n;
    double m_mut_fraction;
    mutable detail::random_engine_type m_e;
    unsigned m_seed;
    unsigned m_verbosity;
    mutable log_type m_log;
    pagmo::bfe m_bfe;
};
} // namespace dcgp
#endif
