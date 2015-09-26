
#include <iostream>
#include <sstream>
#include <random>
#include <limits>
#include <cmath>
#include <stdexcept>

#include "expression.h"
#include "std_overloads.h"


namespace dcgp {

/// Constructor
/** Constructs a d-cgp expression
 *
 * \param[in] n number of inputs (independent variables)
 * \param[in] m number of outputs (dependent variables)
 * \param[in] r number of rows of the cartesian cgp
 * \param[in] c number of columns of the cartesian cgp
 * \param[in] l number of levels-back allowed for the cartesian cgp
 * \param[in] f function set. An std::vector of dcgp::basis_function
 * \param[in] seed seed for the random number generator (initial expression  and mutations depend on this)
 */
expression::expression(unsigned int n,              // n. inputs
                   unsigned int m,                  // n. outputs
                   unsigned int r,                  // n. rows
                   unsigned int c,                  // n. columns
                   unsigned int l,                  // n. levels-back
                   std::vector<basis_function> f,   // functions
                   unsigned int seed                // seed for the pseudo-random numbers
                   ) : m_n(n), m_m(m), m_r(r), m_c(c), m_l(l), m_f(f), m_lb((3 * m_r * m_c) + m_m, 0), m_ub((3 * m_r * m_c) + m_m, 0), m_x((3 * m_r * m_c) + m_m, 0), m_e(seed)
{

    if (n == 0) throw std::invalid_argument("Number of inputs is 0");
    if (m == 0) throw std::invalid_argument("Number of outputs is 0");
    if (c == 0) throw std::invalid_argument("Number of columns is 0");
    if (r == 0) throw std::invalid_argument("Number of rows is 0");
    if (l == 0) throw std::invalid_argument("Number of level-backs is 0");
    if (f.size()==0) throw std::invalid_argument("Number of basis functions is 0");

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

    // We generate a random expression
    for (auto i = 0u; i < m_x.size(); ++i)
    {
        m_x[i] = std::uniform_int_distribution<unsigned int>(m_lb[i], m_ub[i])(m_e);
    }
    update_active();
}

/// Sets the chromosome
/** Sets a new chromosome as genotype for the expression and updates the active nodes and active genes information
 *
 * \param[in] x The new cromosome
 *
 * @throw dcgp::input_error if the chromosome is incompatible with the expression (n.inputs, n.outputs, levels-back, etc.)
 */
void expression::set(const std::vector<unsigned int>& x)
{
    if(!is_valid(x))
    {
        throw std::invalid_argument("Chromosome is incompatible");
    }
    m_x = x;
    update_active();
}


inline unsigned int factorial(unsigned int n)
{
    if (n==0) return 1;
    unsigned int ret = 1;
    for(auto i = 1u; i <= n; ++i)
        ret *= i;
    return ret;
}


/// Computes the derivatives of the expression
/** 
 * Using automated differentiation rules this method returns the derivatives up to a certain order, with respect
 * to one input variable at a given point.
 *
 * \param[in] wrt index of the derivation variable (0,1 ..., m_n)
 * \param[in] order the derivative order we want to compute
 * \param[in] in std::vector containing the point coordinates we want the derivatives be computed at
 *
 * @returns std::vector<std::vector<double> > containing the value of f,f',f'' ..., at the point in
 *
 * @throw std::invalid_argument
 */
std::vector<std::vector<double> > expression::differentiate(unsigned int wrt, unsigned int order, const std::vector<double>& in) const
{  
    if(in.size() != m_n)
    {
        throw std::invalid_argument("Input size is incompatible");
    }
    if(wrt >= m_n)
    {
        throw std::invalid_argument("Derivative id is larger than the independent variable number");
    }
//for (auto i : m_active_nodes) std::cout << " " << i; std::cout << std::endl;
    std::vector<double> dumb(m_m);
    std::vector<std::vector<double> > retval(order+1,dumb);
    std::map<unsigned int, std::vector<double> > node_jet;
    for (auto j =0u; j<=order; ++j)
    {
        for (auto i : m_active_nodes)
        {
            if (i < m_n) 
            {
                //if (j==0) node_jet[i] = std::vector<double>({in[i]});
                if (j==0) node_jet[i].push_back(in[i]);
                else if (j==1) node_jet[i].push_back((i==wrt) ? 1. : 0.);
                else node_jet[i].push_back(0.);
            } else {
                unsigned int idx = (i - m_n) * 3;
                if (j==0) node_jet[i] = std::vector<double>({m_f[m_x[idx]].m_df(node_jet[m_x[idx + 1]], node_jet[m_x[idx + 2]])});
                else node_jet[i].push_back(m_f[m_x[idx]].m_df(node_jet[m_x[idx + 1]], node_jet[m_x[idx + 2]]));
            }
        }
    }
//std::cout << i << ", " << node_jet[i] << std::endl;

    for (auto j = 0u; j<=order; ++j) {
        for (auto i = 0u; i<m_m; ++i)
        {
            retval[j][i] = node_jet[m_x[(m_r * m_c) * 3 + i]][j] * factorial(j);
        }
    }
    return retval;
}

/// Mutates one of the active genes
/** 
 * Mutates exactly one of the active genes
 */
void expression::mutate_active()
{
    unsigned int idx = std::uniform_int_distribution<unsigned int>(0, m_active_genes.size() - 1)(m_e);
    idx = m_active_genes[idx];

    if (m_lb[idx]<m_ub[idx]) // if only one value is allowed for the gene, then we will not do anything as mutation does not apply
    {
        unsigned int new_value = UINT_MAX;
        do 
        {
            new_value = std::uniform_int_distribution<unsigned int>(m_lb[idx], m_ub[idx])(m_e);
        } while (new_value == m_x[idx]);
        m_x[idx] = new_value;
        update_active();
    }
}

/// Validity of a chromosome
/** 
 * Checks if a chromosome (i.e. a sequence of integers) is a valid expression
 * by checking its length and the bounds
 *
 * \param[in] x chromosome 
 */
bool expression::is_valid(const std::vector<unsigned int>& x) const
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


/// Updates the m_active_genes and m_active_nodes data member
void expression::update_active()
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
        // We remove duplicates to avoid processng them and thus having a 2^N complexity
        std::sort( next.begin(), next.end() );
        next.erase( std::unique( next.begin(), next.end() ), next.end() );
        current = next;
        next.clear();
    } while (current.size() > 0);

    // We remove duplicates and keep m_active_nodes sorted
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
std::string expression::human_readable() const
{
    std::ostringstream s;
    s << "d-CGP Expression:\n";
    s << "\tNumber of inputs:\t\t" << m_n << '\n';
    s << "\tNumber of outputs:\t\t" << m_m << '\n';
    s << "\tNumber of rows:\t\t\t" << m_r << '\n';
    s << "\tNumber of columns:\t\t" << m_c << '\n';
    s << "\tNumber of levels-back allowed:\t" << m_l << '\n';
    s << "\n\tResulting lower bounds:\t" << m_lb;
    s << "\n\tResulting upper bounds:\t" << m_ub << '\n';
    s << "\n\tCurrent expression (encoded):\t" << m_x << '\n';
    s << "\tActive nodes:\t\t\t" << m_active_nodes << '\n';
    s << "\tActive genes:\t\t\t" << m_active_genes << '\n';
    s << "\n\tFunction set:\t\t\t" << m_f << '\n';
    return s.str();
}

/// Overload stream operator for dcgp::expression
/**
 * Equivalent to printing expression::human_readable() to stream.
 *
 * @param[out] s std::ostream to which the problem will be streamed.
 * @param[in] p dcgp::expression to be inserted into the stream.
 *
 * @return reference to s.
 */
std::ostream &operator<<(std::ostream &s, const expression &en)
{
    s << en.human_readable();
    return s;
}

} // end of namespace dcgp

