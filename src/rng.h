#include <random>
namespace dcgp {
class rng {
public:
	static unsigned int get_seed();
private:
	static std::random_device m_rdev;
};
}