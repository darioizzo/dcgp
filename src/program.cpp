
#include <iostream>
#include <sstream>
#include <random>
#include <limits>

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
                   double tol,                      // tolerance for the fitness (only when hits based)
                   unsigned int seed                // seed for the pseudo-random numbers
                   ) : m_n(n), m_m(m), m_r(r), m_c(c), m_l(l), m_f(f), m_lb((3 * m_r * m_c) + m_m, 0), m_ub((3 * m_r * m_c) + m_m, 0), m_x((3 * m_r * m_c) + m_m, 0), m_tol(tol), m_e(seed)
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
std::cout << "I AM HERE" << std::endl;
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

/// Gets the program
const std::vector<unsigned int>&  program::get() const
{
    return m_x;
}

/// Gets the active genes
const std::vector<unsigned int>&  program::get_active_genes() const
{
    return m_active_genes;
}

/// Gets the active nodes
const std::vector<unsigned int>&  program::get_active_nodes() const
{
    return m_active_nodes;
}

/// Computes the error of the program in approximating some given data
double program::fitness(const std::vector<std::vector<double> >& in_des, const std::vector<std::vector<double> >& out_des, fitness_type type) const
{
    double retval = 0.;
    std::vector<double> out_real;

    if (in_des.size() != out_des.size())
    {
        throw input_error("Size of the input vector must be the size of the output vector");
    }

    for (auto i = 0u; i < in_des.size(); ++i)
    {
        out_real = compute_f(in_des[i]);
        if (type == fitness_type::ERROR_BASED)
        {
            for (auto j = 0u; j < out_real.size(); ++j)
            {
                retval += 1.0 / (1.0 + fabs(out_des[i][j] - out_real[i]));
            }
        } else if (type == fitness_type::HITS_BASED){
            for (auto j = 0u; j < out_real.size(); ++j)
            {
                if (fabs(out_des[i][j] - out_real[i]) < m_tol) retval += 1.0;
            }
        }
    }
    return retval;
}

void program::mutate()
{
    unsigned int idx = std::uniform_int_distribution<unsigned int>(0, m_active_genes.size() - 1)(m_e);
    idx = m_active_genes[idx];
    unsigned int new_value = UINT_MAX;
    do 
    {
        new_value = std::uniform_int_distribution<unsigned int>(m_lb[idx], m_ub[idx])(m_e);

    } while (new_value == m_x[idx]);
    m_x[idx] = new_value;
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
    std::vector<unsigned int> current(m_m), next;
    m_active_nodes.clear();
    // At the beginning current contains only the output nodes connections
    for (auto i = 0u; i < m_m; ++i) {
        current[i] = m_x[3 * m_r * m_c + i];
    }

    do
    {
        m_active_nodes.insert(m_active_nodes.end(), current.begin(), current.end());
std::cout << "D: " << current.size() << std::endl;
        for (auto node_id : current)
        {
            if (node_id >=m_n) // we insert the input nodes connections as they do not have any
            {
                next.push_back(m_x[(node_id - m_n) * 3 + 1]);
                next.push_back(m_x[(node_id - m_n) * 3 + 2]);
            }
            else{
                m_active_nodes.push_back(node_id);
            }
        }
        current = next;
        next.clear();
    } while (current.size() > 0);
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

