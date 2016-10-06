#include <vector>
#include <random>
#include <iostream>

#include "../../include/expression.hpp"

struct es_params {
    unsigned int m_childs;
    std::string m_mutation_type;
    unsigned int m_n;
    double m_mut_prob;
    unsigned int m_gen;
};

// Evolves the expression ex to fit the in-out supervised data
void es(const std::vector<std::vector<double> >& in, const std::vector<std::vector<double> >& out, dcgp::expression<double> &ex, const es_params& p) {
    // Random seed
    std::random_device rd;
    std::default_random_engine re(rd());

    double best_fit = 1e32;
    std::vector<double> newfits(p.m_childs,0.);
    std::vector<std::vector<unsigned int> > newchromosomes(p.m_childs);
    std::vector<unsigned int> best_chromosome(ex.get());
    unsigned int gen = 0;

    do
    {
        gen++;
        for (auto i = 0u; i < newfits.size(); ++i) {
            ex.set(best_chromosome);
            if (p.m_mutation_type == "active") {
                ex.mutate_active(p.m_n);
            } else {
                std::vector<unsigned int> tbm;
                for (auto i = 0u; i < best_chromosome.size(); ++ i) {
                    if (std::uniform_real_distribution<double>(0, 1)(re) < p.m_mut_prob) tbm.push_back(i);
                }
                ex.mutate(tbm);
            }
            newfits[i] = quadratic_error(ex, in, out);
            newchromosomes[i] = ex.get();
        }

        for (auto i = 0u; i < newfits.size(); ++i) {
            if (newfits[i] <= best_fit) {
                if (newfits[i] != best_fit) {
                    std::cout << "New best found: gen: " << std::setw(7) << gen << "\t value: " << newfits[i] << std::endl;
                }
                best_fit = newfits[i];
                best_chromosome = newchromosomes[i];
                ex.set(best_chromosome);
            }
        }
    } while (best_fit > 1e-3 && gen < p.m_gen);
    ex.set(best_chromosome);
    std::cout << "Number of generations: " << gen << std::endl;
}
