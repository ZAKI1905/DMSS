#ifndef Zaki_Math_IntegralTable_H
#define Zaki_Math_IntegralTable_H

#include "Zaki/Math/Math_Core.hpp"

//==============================================================
namespace Zaki::Math
{
//==============================================================
/// I_1 integral with integrand :
///     x^2 / sqrt[ x^2 + a^2 ]
double I_1(const Range<double>& limits, const std::vector<double>& pars) ;

/// I_2 integral with integrand :
///     k^2 * sqrt[ k^2 + a^2 ]
double I_2(const Range<double>& limits, const std::vector<double>& pars) ;

/// I_3 integral with integrand :
///     k^4 / sqrt[ k^2 + a^2 ]
double I_3(const Range<double>& limits, const std::vector<double>& pars) ;

//==============================================================
//--------------------------------------------------------------
} // End of namespace Zaki::Math
//--------------------------------------------------------------

#endif /*Zaki_Math_IntegralTable_H*/