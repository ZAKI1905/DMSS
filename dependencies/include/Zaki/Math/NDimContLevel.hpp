#ifndef Zaki_Math_NDimContLevel_H
#define Zaki_Math_NDimContLevel_H

// Zaki::Math
#include "Zaki/Math/Math_Core.hpp"

//==============================================================
namespace Zaki::Math
{

//==============================================================
class NDimContLevel
{
  private:

    double PSeries(const double n) const;

    double P(const double n) const ;

    double tolerance = 0.00001 ;
    int dim = 1 ;
    double conf_level = 1 ; // 1 sigma


  public:

    /// Default constructor
    NDimContLevel() ;

    /// First constructor
    NDimContLevel(const int in_Dim, const double in_CL) 
     : dim(in_Dim), conf_level(in_CL) {} ;

    /// Sets the tolerance for the solver
    void SetTolerance(const double) ;
    
    /// Sets the dimension
    void SetDim(const int) ;
    
    /// Sets the Confidence Level C.L.
    void SetCL(const double) ;

    struct Solution
    {
      int Dim ; double CL, P, Err;
      double Val, Up;
    };

    /// Solves for Delta(L) corresponding to the confidence level
    Solution Solve() ;

};

//==============================================================
//                  ostream << overloading
std::ostream& operator << (std::ostream &output,
const Zaki::Math::NDimContLevel::Solution&) ;
//==============================================================

//--------------------------------------------------------------
} // End of namespace Zaki::Math
//--------------------------------------------------------------

#endif /*Zaki_Math_NDimContLevel_H*/