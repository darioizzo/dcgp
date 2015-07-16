#include <boost/lexical_cast.hpp>

#define DCGP_MAX_OUTPUT_LENGTH 20
namespace std
{
	/// Overload stream insertion operator for std::vector<T>. It will only output the first
	/// DCGP_MAX_OUTPUT_LENGTH elements.
	template < class T >
	inline ostream &operator<<(ostream &os, const vector<T> &v)
	{
		typename vector<T>::size_type len = v.size();
		if (len < DCGP_MAX_OUTPUT_LENGTH) 
		{
			os << '[';
			for (typename std::vector<T>::size_type i = 0; i < v.size(); ++i) {
				os << boost::lexical_cast<std::string>(v[i]);
				if (i != v.size() - 1) {
					os << ", ";
				}
			}
			os << ']';
		} else {
			os << '[';
			for (typename std::vector<T>::size_type i = 0; i < DCGP_MAX_OUTPUT_LENGTH; ++i) {
				os << boost::lexical_cast<std::string>(v[i]) << ", ";
			}
			os << " ... ]";
		}
		return os;
	}

	
}
#undef DCGP_MAX_OUTPUT_LENGTH