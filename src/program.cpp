
#include <iostream>
#include <sstream>
#include <random>

#include "program.h"
#include "std_overloads.h"


namespace dcgp {

/// Constructor
/** Constructs a d-cgp program
 *
 * \param[in] n number of inputs (independent variables)
 * \param[in] m number of outputs (dependent variables)
 * \param[in] r number of rows of the cartesian cgp
 * \param[in] c number of columns of the cartesian cgp
 * \param[in] l number of levels-back allowed for the cartesian cgp
 * \param[in] f functions allowed. They are of type dcgp::basis_function
 */
program::program(unsigned int n,                    // n. inputs
                   unsigned int m,                  // n. outputs
                   unsigned int r,                  // n. rows
                   unsigned int c,                  // n. columns
                   unsigned int l,                  // n. levels-back
                   std::vector<basis_function> f,   // functions
                   unsigned int seed                // seed for the pseudo-random numbers
                   ) : m_n(n), m_m(m), m_r(r), m_c(c), m_l(l), m_f(f), m_lb((3 * m_r * m_c) + m_m, 0), m_ub((3 * m_r * m_c) + m_m, 0), m_x((3 * m_r * m_c) + m_m, 0), m_e(seed)
{

    if (n == 0) throw input_error("Number of inputs is 0");
    if (m == 0) throw input_error("Number of outputs is 0");
    if (c == 0) throw input_error("Number of columns is 0");
    if (r == 0) throw input_error("Number of rows is 0");
    if (l == 0) throw input_error("Number of level-backs is 0");
    if (f.size()==0) throw input_error("Number of basis functions is 0");

    // Bounds for the function genes
    for (auto i = 0u; i < (3 * m_r * m_c); i+=3) {
        m_ub[i] = f.size() - 1;
    }

    // Bounds for the output genes
    for (auto i = 3u * m_r * m_c; i < m_ub.size(); ++i) {
        m_ub[i] = m_n + m_r * m_c - 1;
        if (m_l <= m_c) {
            m_lb[i] = m_n + m_r * (m_c - m_l);
        }
    }

    // Bounds for the node connection genes 
    for (auto i = 0u; i < m_c; ++i) {
        for (auto j = 0u; j < m_r; ++j) {
            m_ub[((i * m_r) + j) * 3 + 1] = m_n + i * m_r - 1;
            m_ub[((i * m_r) + j) * 3 + 2] = m_n + i * m_r - 1;
            if (i >= m_l) {
                m_lb[((i * m_r) + j) * 3 + 1] = m_n + m_r * (i - m_l);
                m_lb[((i * m_r) + j) * 3 + 2] = m_n + m_r * (i - m_l);
            }
        }
    }

    // We generate a random program
    for (auto i = 0u; i < m_x.size(); ++i)
    {
        m_x[i] = std::uniform_int_distribution<unsigned int>(m_lb[i], m_ub[i])(m_e);
    }
    update_active();
}

/// Sets the pregram
void program::set(const std::vector<unsigned int>& x)
{
    if(!is_valid(x))
    {
        throw input_error("Chromosome is incompatible");
    }
    m_x = x;
    update_active();
}

/// Validity of a chromosome
/** Checks if a chromosome (i.e. a sequence of integers) is a valid program
 *
 * \param[in] x chromosome 
 */
bool program::is_valid(const std::vector<unsigned int>& x) const
{
    // Checking for length
    if (x.size() != m_lb.size()) {
        return false;
    }

    // Checking for bounds on all cenes
    for (auto i = 0u; i < x.size(); ++i) {
        if ((x[i] > m_ub[i]) || (x[i] < m_lb[i])) {
            return false;
        }
    }
    return true;
}


/// Computes which nodes actually need evaluation
void program::update_active()
{
    assert(m_x.size() == m_lb.size());

    // First we update the active nodes
    m_active_nodes.resize((m_c * m_r) * 2 + m_m);

    // We start with the output nodes and their dependencies
    for (auto i = 0u; i < m_m; ++i) {
        m_active_nodes[(m_c * m_r) * 2 + i] = m_x[(m_c * m_r) * 3 + i];
    }

    // We start with the output nodes and their dependencies
    for (auto i = 0u; i < m_c * m_r; ++i) {
        m_active_nodes[2*i] = m_x[3 * i + 1];
        m_active_nodes[2*i + 1] = m_x[3 * i + 2];
    }
    std::sort( m_active_nodes.begin(), m_active_nodes.end() );
    m_active_nodes.erase( std::unique( m_active_nodes.begin(), m_active_nodes.end() ), m_active_nodes.end() );

    // Then the active genes
    m_active_genes.clear();
    for (auto i = 0u; i<m_active_nodes.size(); ++i) 
    {
        if (m_active_nodes[i] >= m_n) 
        {
            unsigned int idx = (m_active_nodes[i] - m_n) * 3;
            m_active_genes.push_back(idx);
            m_active_genes.push_back(idx + 1);
            m_active_genes.push_back(idx + 2);
        }
    }
    for (auto i = 0u; i<m_m; ++i) 
    {
        m_active_genes.push_back(m_r * m_c * 3 + i);
    }
}

/// Return human readable representation of the problem.
/**
 * Will return a formatted string containing a human readable representation of the class
 *
 * @return std::string containing a human-readable representation of the problem.
 */
std::string program::human_readable() const
{
    std::ostringstream s;
    s << "CGP Encoding:\n";
    s << "\tNumber of inputs:\t\t" << m_n << '\n';
    s << "\tNumber of outputs:\t\t" << m_m << '\n';
    s << "\tNumber of rows:\t\t\t" << m_r << '\n';
    s << "\tNumber of columns:\t\t" << m_c << '\n';
    s << "\tNumber of levels-back allowed:\t" << m_l << '\n';
    s << "\n\tResulting lower bounds:\t" << m_lb;
    s << "\n\tResulting upper bounds:\t" << m_ub << '\n';
    s << "\n\tCurrent program (encoded):\t" << m_x << '\n';
    s << "\n\tActive nodes:\t" << m_active_nodes << '\n';
    s << "\n\tActive genes:\t" << m_active_genes << '\n';
    return s.str();
}

/// Overload stream operator for problem::base.
/**
 * Equivalent to printing program::human_readable() to stream.
 *
 * @param[out] s std::ostream to which the problem will be streamed.
 * @param[in] p dcgp::program to be inserted into the stream.
 *
 * @return reference to s.
 */
std::ostream &operator<<(std::ostream &s, const program &en)
{
    s << en.human_readable();
    return s;
}

} // end of namespace dcgp

