#include "rng.h"

namespace dcgp {

std::random_device rng::m_rdev;

unsigned int rng::get_seed() {return m_rdev();}
}