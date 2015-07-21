#include <vector>
#include "expression.h"

namespace dcgp {
    enum fitness_type { 
        ERROR_BASED,    // fitness is sum_ij [1 / (1 + err_ij)] 
        HITS_BASED      // fitness is the number of components across the output data which are within a tolerance
        };   

    /// Computes the error of the expression in approximating some given data
    double simple_data_fit(const dcgp::expression& ex, 
        const std::vector<std::vector<double> >& in_des, 
        const std::vector<std::vector<double> >& out_des, 
        fitness_type type = fitness_type::ERROR_BASED,
        double tol = 1e-10);
}
