#include <boost/python.hpp>
#include <pygmo/algorithm_exposition_suite.hpp>
#include <pygmo/register_ap.hpp>
#include <sstream>
#include <string>
#include <vector>

#include <dcgp/gym.hpp>

// See: https://docs.scipy.org/doc/numpy/reference/c-api.array.html#importing-the-api
// In every cpp file We need to make sure this is included before everything else,
// with the correct #defines.
#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL dcgpy_ARRAY_API
#include "common_utils.hpp"

#include "docstrings.hpp"

namespace bp = boost::python;
using namespace dcgp;

using gym_ptr = void (*)(std::vector<std::vector<double>> &, std::vector<std::vector<double>> &);
template <gym_ptr F>
inline void expose_data_from_the_gym(const std::string &name, const std::string &docstring = std::string{"ds"})
{
    bp::def(
        name.c_str(),
        +[]() {
            std::vector<std::vector<double>> points, labels;
            F(points, labels);
            return bp::make_tuple(dcgpy::vvector_to_ndarr<double>(points), dcgpy::vvector_to_ndarr<double>(labels));
        },
        docstring.c_str());
}

namespace dcgpy
{

void expose_symbolic_regression()
{
    // Making data from the gym available in python
    expose_data_from_the_gym<&gym::generate_koza_quintic>("generate_koza_quintic", generate_koza_quintic_doc());
    // From Our paper
    expose_data_from_the_gym<&gym::generate_P1>("generate_P1", generate_P1_doc());
    expose_data_from_the_gym<&gym::generate_P2>("generate_P2", generate_P2_doc());
    expose_data_from_the_gym<&gym::generate_P3>("generate_P3", generate_P3_doc());
    expose_data_from_the_gym<&gym::generate_P4>("generate_P4", generate_P4_doc());
    expose_data_from_the_gym<&gym::generate_P5>("generate_P5", generate_P5_doc());
    expose_data_from_the_gym<&gym::generate_P6>("generate_P6", generate_P6_doc());
    expose_data_from_the_gym<&gym::generate_P7>("generate_P7", generate_P7_doc());
    // From Vladi paper
    expose_data_from_the_gym<&gym::generate_kotanchek>("generate_kotanchek", generate_kotanchek_doc());
    expose_data_from_the_gym<&gym::generate_salutowicz>("generate_salutowicz", generate_salutowicz_doc());
    expose_data_from_the_gym<&gym::generate_salutowicz2d>("generate_salutowicz2d", generate_salutowicz2d_doc());
    expose_data_from_the_gym<&gym::generate_uball5d>("generate_uball5d", generate_uball5d_doc());
    expose_data_from_the_gym<&gym::generate_ratpol3d>("generate_ratpol3d", generate_ratpol3d_doc());
    expose_data_from_the_gym<&gym::generate_sinecosine>("generate_sinecosine", generate_sinecosine_doc());
    expose_data_from_the_gym<&gym::generate_ripple>("generate_ripple", generate_ripple_doc());
    expose_data_from_the_gym<&gym::generate_ratpol2d>("generate_ratpol2d", generate_ratpol2d_doc());
    // NIST data
    expose_data_from_the_gym<&gym::generate_chwirut1>("generate_chwirut1", generate_chwirut1_doc());
    expose_data_from_the_gym<&gym::generate_chwirut2>("generate_chwirut2", generate_chwirut2_doc());
    expose_data_from_the_gym<&gym::generate_daniel_wood>("generate_daniel_wood", generate_daniel_wood_doc());
    expose_data_from_the_gym<&gym::generate_gauss1>("generate_gauss1", generate_gauss1_doc());
    expose_data_from_the_gym<&gym::generate_kirby2>("generate_kirby2", generate_kirby2_doc());
    expose_data_from_the_gym<&gym::generate_lanczos2>("generate_lanczos2", generate_lanczos2_doc());
    expose_data_from_the_gym<&gym::generate_misra1b>("generate_misra1b", generate_misra1b_doc());
}
} // namespace dcgpy