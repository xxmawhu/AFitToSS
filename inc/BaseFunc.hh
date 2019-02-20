#ifndef Roo_Propagator_HH
#define Roo_Propagator_HH

#include <iostream>
#include <math.h>

#if defined(USEROOT) || defined(__CINT__)
#else
#endif
#include "TComplex.h"
#include <vector>
using std::vector;
namespace Func{
Double_t dot(const Double_t p4_1[4],
        const Double_t p4_2[4]);
Double_t det(  
        const Double_t a[4],
        const Double_t b[4],
        const Double_t c[4],
        const Double_t d[4]
        );
};

#endif
