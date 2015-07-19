#ifndef DCGP_EXCEPTIONS_H
#define DCGP_EXCEPTIONS_H

#include <string>

namespace dcgp {

/// An exception originating from derivative computations.
class derivative_error : public std::exception {
public:
    derivative_error(std::string message) : m_message(message) {}
    ~derivative_error() throw() {}
    const char* what() const throw() {  
        return m_message.c_str();
    }
private:
    std::string m_message;
};

/// An exception originated by detected inconsistent inputs.
class input_error : public std::exception {
public:
    input_error(std::string message) : m_message(message) {}
    ~input_error() throw() {}
    const char* what() const throw() {  
        return m_message.c_str();
    }
private:
    std::string m_message;
};

} // end of namespace dcgp

#endif // DCGP_DCGP_H
