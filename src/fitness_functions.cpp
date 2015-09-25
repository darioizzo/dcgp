#include <cmath>
#include <vector>
#include <stdexcept>

#include "fitness_functions.h"
#include "expression.h"


namespace dcgp {
    /// Computes the error of the expression in approximating some given data
    double simple_data_fit(const expression& ex, 
        const std::vector<std::vector<double> >& in_des, 
        const std::vector<std::vector<double> >& out_des, 
        fitness_type type,
        double tol) 
    {
        double retval = 0.;
        std::vector<double> out_real;

        if (in_des.size() != out_des.size())
        {
            throw std::invalid_argument("Size of the input vector must be the size of the output vector");
        }

        for (auto i = 0u; i < in_des.size(); ++i)
        {
            out_real = ex(in_des[i]);
            if (type == fitness_type::ERROR_BASED)
            {
                for (auto j = 0u; j < out_real.size(); ++j)
                {
                    if (std::isfinite(out_real[j]))
                    {
                        retval += 1.0 / (1.0 + fabs(out_des[i][j] - out_real[j]));
                    }
                }
            } else if (type == fitness_type::HITS_BASED){
                for (auto j = 0u; j < out_real.size(); ++j)
                {
                    if (std::isfinite(out_real[j]))
                    {
                        if (fabs(out_des[i][j] - out_real[j]) < tol) retval += 1.0;
                    }
                }
            }
        }

        return retval;
    }
}
