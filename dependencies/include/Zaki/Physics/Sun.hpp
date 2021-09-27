//
//  Sun.hpp
//  zaki
//
//  Created by Mohammadreza Zakeri on 5/2/20.
//

#ifndef Zaki_Physics_Sun_hpp
#define Zaki_Physics_Sun_hpp

#include "Zaki/Physics/Coordinate.hpp"

//--------------------------------------------------------------
namespace Zaki::Physics
{
    
//==============================================================
//class Sun
//{
    
//public:
    /// Sun position in GEI coordinates
    Zaki::Physics::GEICoord GetSunPos(const double&) ;

    /// Earth position in Heliocentric Ecliptic Coordinate
    Zaki::Physics::HEclipticCoord GetEarthPos(const double&) ;
    
    ///  The ecliptic longitude of the Sun
    double GetSunEclipticLong(const double&) ;
//};
//==============================================================

} // End of namespace Zaki::Physics
//--------------------------------------------------------------
#endif /* Zaki_Physics_Sun_hpp */
