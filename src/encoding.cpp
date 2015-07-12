
#include <iostream>

#include "encoding.h"


namespace dcgp {

/// Constructor
/** Constructs a differential genetic programming encoding
 *
 * \param[in] n number of inputs (independent variables)
 * \param[in] m number of outputs (dependent variables)
 * \param[in] c number of columns of the cartesian cgp
 * \param[in] r number of rows of the cartesian cgp
 * \param[in] l number of levels-back allowed for the cartesian cgp
 * \param[in] f functions allowed. They are of type dcgp::basis_function
 */
encoding::encoding(unsigned int n,                  // n. inputs
                   unsigned int m,                  // n. outputs
                   unsigned int c,                  // n. columns
                   unsigned int r,                  // n. rows
                   unsigned int l,                  // n. levels-back
                   std::vector<basis_function> f    // functions
                   ) : m_n(n), m_m(m), m_c(c), m_r(r), m_l(l), m_f(f), m_lb((3 * m_r * m_c) + m_m, 0), m_ub((3 * m_r * m_c) + m_m, 0)
{
    // Bounds for the function genes
    for (auto i = 0; i < (3 * m_r * m_c); i+=3) {
        m_ub[i] = f.size() - 1;
    }

    // Bounds for the output genes
    for (auto i = 3 * m_r * m_c; i < m_ub.size(); ++i) {
        m_ub[i] = m_n + m_r * m_c - 1;
        if (m_l <= m_c) {
            m_lb[i] = m_n + m_r * (m_c - m_l);
        }
    }

    // Bounds for the node connection genes 
    for (auto i = 0; i < m_c; ++i) {
        for (auto j = 0; j < m_r; ++j) {
            m_ub[((i * m_r) + j) * 3 + 1] = m_n + i * m_r - 1;
            m_ub[((i * m_r) + j) * 3 + 2] = m_n + i * m_r - 1;
            if (i >= m_l) {
                m_lb[((i * m_r) + j) * 3 + 1] = m_n + m_r * (i - m_l);
                m_lb[((i * m_r) + j) * 3 + 2] = m_n + m_r * (i - m_l);
            }
        }
    }

}

/// Validity of a chromosome
/** Checks if a chromosome (i.e. a sequence of integers) can be decoded with this encoding
 *
 * \param[in] x chromosome 
 * \param[in] verbose if true prints additinal information in case the chromosome is not valid
 */
bool encoding::is_valid(std::vector<unsigned int> x, bool verbose)
{
    // Checking for length
    if (x.size() != m_lb.size()) {
        if (verbose) std::cout << "Chromosome dimension is wrong" << std::endl;
        return false;
    }

    for (auto i = 0u; i < x.size(); ++i) {
        if ((x[i] > m_ub[i]) || (x[i] < m_lb[i])) {
            if (verbose) {
                std::cout << "Out of bounds: \n";
                std::cout << "lb = [";
                for (auto i = 0u; i<m_lb.size(); ++i) {
                    std::cout << m_lb[i] << " ";
                }
                std::cout << "]\nub = [";
                for (auto i = 0u; i<m_ub.size(); ++i) {
                    std::cout << m_ub[i] << " ";
                }
                std::cout << "]\nx = [";
                for (auto i = 0u; i<x.size(); ++i) {
                    std::cout << x[i] << " ";
                }
                std::cout << "]\n";
            }
            return false;
        }
    }
    return true;
}

} // end of namespace dcgp

