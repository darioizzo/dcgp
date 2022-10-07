#ifndef DCGP_SYMENGINE_H
#define DCGP_SYMENGINE_H

// patch to make this compile in clang-cl
#if defined(_MSC_VER) && defined(__clang__)
#define and &&
#define or ||
#define not !
#endif

#include <symengine/expression.h>

#if defined(_MSC_VER) && defined(__clang__)
#undef and
#undef or
#undef not
#endif

// patch to avoid flint defining access as  _access in windows and messing with boost
#if defined(access)
#undef access
#endif

#endif
