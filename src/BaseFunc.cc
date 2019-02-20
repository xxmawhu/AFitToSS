#include "TComplex.h"
#include "math.h"
#include "BaseFunc.hh"
Double_t Func::dot(const Double_t p4_1[4],
        const Double_t p4_2[4])
{
    return p4_1[0]*p4_2[0] -p4_1[1]*p4_2[1] - p4_1[2]*p4_2[2] 
        -p4_1[3]*p4_2[3] ;
}
Double_t Func::det(  
        const Double_t a[4],
        const Double_t b[4],
        const Double_t c[4],
        const Double_t d[4]
        )
{
    return 
        a[1]*b[3]*c[2]*d[0] - a[1]*b[2]*c[3]*d[0] - a[0]*b[3]*c[2]*d[1] + 

        a[0]*b[2]*c[3]*d[1] - a[1]*b[3]*c[0]*d[2] + a[0]*b[3]*c[1]*d[2] + 

        a[1]*b[0]*c[3]*d[2] - a[0]*b[1]*c[3]*d[2] + 

        a[3]*(-(b[1]*c[2]*d[0]) + b[0]*c[2]*d[1] + b[2]*(c[1]*d[0] - c[0]*d[1]) + 

                b[1]*c[0]*d[2] - b[0]*c[1]*d[2]) + a[1]*b[2]*c[0]*d[3] - 

        a[0]*b[2]*c[1]*d[3] - a[1]*b[0]*c[2]*d[3] + a[0]*b[1]*c[2]*d[3] + 

        a[2]*(-(b[3]*c[1]*d[0]) + b[1]*c[3]*d[0] + b[3]*c[0]*d[1] - 

     b[0]*c[3]*d[1] - b[1]*c[0]*d[3] + b[0]*c[1]*d[3]) 
      ;
}
