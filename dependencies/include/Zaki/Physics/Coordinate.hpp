#ifndef Zaki_Physics_Coordinate_H
#define Zaki_Physics_Coordinate_H

// Root
#include <Math/Vector3D.h>

#include "Zaki/Physics/DateTime.hpp"
#include "Zaki/Physics/Constants.hpp"

//--------------------------------------------------------------
namespace Zaki::Physics
{

//==============================================================
struct Coord3D
{
    double x,y,z ;

    Coord3D() {}
    Coord3D(const double in_x, const double in_y, const double in_z) 
      : x(in_x), y(in_y), z(in_z) {}

    bool operator==(const Coord3D& rhs) const
    {
        return ( x == rhs.x && y == rhs.y && z == rhs.z ) ;
    }

    bool operator>(const Coord3D& rhs) const
    {
        return ( XYDist2({0,0,0}) > rhs.XYDist2({0,0,0})) ;
    }

    bool operator<(const Coord3D& rhs) const
    {
        return ( XYDist2({0,0,0}) < rhs.XYDist2({0,0,0})) ;
    }

    double XYDist2(const Coord3D& rhs) const
    {
        return (x - rhs.x)*(x - rhs.x)
                + (y - rhs.y)*(y - rhs.y) ;
    }

    std::string Str() const ;

  // coords operator+(const coords& rhs) const
  // {
  //   return {x + rhs.x, y + rhs.y, z + rhs.z } ;
  // }

  // coords operator-(const coords& rhs) const
  // {
  //   return {x - rhs.x, y - rhs.y, z - rhs.z } ;
  // }
  
  // coords operator*(const double& fac) const
  // {
  //   return {fac*x, fac*y, fac*z } ;
  // }
};

//==============================================================
// Yaw, Pitch, Roll Coordinates
struct YPR
{
  double Yaw, Pitch, Roll ;

  YPR Rad() const
  {
    return {Yaw*DEG_2_RAD, Pitch*DEG_2_RAD, Roll*DEG_2_RAD} ;
  }
};

    class GEclipticCoord ;
    class HEclipticCoord ;
    class GEICoord ;
//==============================================================
class CoordBase
{
//    friend std::ostream& operator << (std::ostream &, const CoordBase&) ;
    
protected:
    // Rectangular
    ROOT::Math::XYZVector xyz ;
    
public:
    
    CoordBase(const ROOT::Math::XYZVector& in_xyz) { xyz = in_xyz; }
    CoordBase(const double& in_x, const double& in_y, const double& in_z)
    { xyz = {in_x, in_y, in_z}; }
    
    virtual ~CoordBase() {} ;
    
    ROOT::Math::XYZVector XYZ() const { return xyz ;}
    // .....................................
    // Operators
//    bool operator==(const CoordBase& other)
//    { return xyz == other.xyz; }
//
//    CoordBase operator+(const CoordBase& other)
//    { return {xyz + other.xyz}; }
//
//    CoordBase operator-(const CoordBase& other)
//    { return {xyz - other.xyz}; }
//
//    CoordBase operator*(const int& other)
//    { return {other*xyz}; }
//
//    CoordBase operator*(const double& other)
//    { return {other*xyz}; }
//
//    CoordBase operator-()
//    { return {-xyz}; }
    // .....................................
    
    /// Convert to Geocentric Ecliptic coordinates
    virtual GEclipticCoord ToGEcliptic(const Zaki::Physics::Date&) = 0;
    
    /// Convert to Geocentric Ecliptic coordinates
    /// overloaded for UnixTDays input
    virtual GEclipticCoord ToGEcliptic(const double&) = 0;
    
    /// Convert to Heliocentric Ecliptic coordinates
    virtual HEclipticCoord ToHEcliptic(const Zaki::Physics::Date&) = 0;
    
    /// Convert to Heliocentric Ecliptic coordinates
    /// overloaded for UnixTDays input
    virtual HEclipticCoord ToHEcliptic(const double&) = 0;
    
    /// Convert to Geocentric Equatorial coordinates
    virtual GEICoord ToGEI(const Zaki::Physics::Date&) = 0;
    
    /// Convert to Geocentric Equatorial coordinates
    /// overloaded for UnixTDays input
    virtual GEICoord ToGEI(const double&) = 0;
    
};
//==============================================================
// Geocentric Ecliptic Coordinate
class GEclipticCoord : public CoordBase
{
    friend std::ostream& operator << (std::ostream &, const GEclipticCoord&) ;
    
public:
    
    GEclipticCoord(const ROOT::Math::XYZVector& in_xyz)
    : CoordBase(in_xyz) {}
    
    GEclipticCoord(const double& in_x, const double& in_y, const double& in_z)
    :  CoordBase(in_x, in_y, in_z) {}
    
    double lambda() const ;
    double beta()   const ;
    double r()      const ;
    
    // .....................................
    // Operators
    bool operator==(const GEclipticCoord& other)
    { return xyz == other.xyz; }
    
    GEclipticCoord operator+(const GEclipticCoord& other)
    { return {xyz + other.xyz}; }
    
    GEclipticCoord operator-(const GEclipticCoord& other)
    { return {xyz - other.xyz}; }
    GEclipticCoord operator*(const int& other)
    { return {other*xyz}; }
    GEclipticCoord operator*(const double& other)
    { return {other*xyz}; }
    GEclipticCoord operator-()
    { return {-xyz}; }
    // .....................................
    
    /// Convert to Geocentric Ecliptic coordinates
    virtual GEclipticCoord ToGEcliptic(const Zaki::Physics::Date&) override;
    
    /// Convert to Geocentric Ecliptic coordinates
    /// overloaded for UnixTDays input
    virtual GEclipticCoord ToGEcliptic(const double&) override ;
    
    /// Convert to Heliocentric Ecliptic coordinates
    virtual HEclipticCoord ToHEcliptic(const Zaki::Physics::Date&) override ;
    
    /// Convert to Heliocentric Ecliptic coordinates
    /// overloaded for UnixTDays input
    virtual HEclipticCoord ToHEcliptic(const double&) override ;
    
    /// Convert to Geocentric Equatorial coordinates
    virtual GEICoord ToGEI(const Zaki::Physics::Date&) override ;
    
    /// Convert to Geocentric Equatorial coordinates
    /// overloaded for UnixTDays input
    virtual GEICoord ToGEI(const double&) override ;
};

//==============================================================
// Heliocentric Ecliptic Coordinate
class HEclipticCoord : public CoordBase
{
    friend std::ostream& operator << (std::ostream &, const HEclipticCoord&) ;
    
public:
    HEclipticCoord(const ROOT::Math::XYZVector& in_xyz)
    : CoordBase(in_xyz) {}
    
    HEclipticCoord(const double& in_x, const double& in_y, const double& in_z)
    :  CoordBase(in_x, in_y, in_z) {}
    
    double l() const ;
    double b() const ;
    double r() const ;
    
    // .....................................
    // Operators
    bool operator==(const HEclipticCoord& other)
    { return xyz == other.xyz; }
    
    HEclipticCoord operator+(const HEclipticCoord& other)
    { return {xyz + other.xyz}; }
    
    HEclipticCoord operator-(const HEclipticCoord& other)
    { return {xyz - other.xyz}; }
    HEclipticCoord operator*(const int& other)
    { return {other*xyz}; }
    HEclipticCoord operator*(const double& other)
    { return {other*xyz}; }
    HEclipticCoord operator-()
    { return {-xyz}; }
    // .....................................
    
    /// Convert to Geocentric Ecliptic coordinates
    virtual GEclipticCoord ToGEcliptic(const Zaki::Physics::Date&) override;
    
    /// Convert to Geocentric Ecliptic coordinates
    /// overloaded for UnixTDays input
    virtual GEclipticCoord ToGEcliptic(const double&) override ;
    
    /// Convert to Heliocentric Ecliptic coordinates
    virtual HEclipticCoord ToHEcliptic(const Zaki::Physics::Date&) override ;
    
    /// Convert to Heliocentric Ecliptic coordinates
    /// overloaded for UnixTDays input
    virtual HEclipticCoord ToHEcliptic(const double&) override ;
    
    /// Convert to Geocentric Equatorial coordinates
    virtual GEICoord ToGEI(const Zaki::Physics::Date&) override ;
    
    /// Convert to Geocentric Equatorial coordinates
    /// overloaded for UnixTDays input
    virtual GEICoord ToGEI(const double&) override ;
};
    
//==============================================================
// Geocentric Equatorial Inertial("GEI") Coordinate
class GEICoord : public CoordBase
{

    friend std::ostream& operator << (std::ostream &, const GEICoord&) ;
    
public:
    
    GEICoord(const ROOT::Math::XYZVector& in_xyz)
        : CoordBase(in_xyz) {}
    
    GEICoord(const double& in_x, const double& in_y, const double& in_z)
    :  CoordBase(in_x, in_y, in_z) {}

    // .....................................
    // Operators
    bool operator==(const GEICoord& other)
    { return xyz == other.xyz; }

    GEICoord operator+(const GEICoord& other)
    { return {xyz + other.xyz}; }

    GEICoord operator-(const GEICoord& other)
    { return {xyz - other.xyz}; }
    GEICoord operator*(const int& other)
    { return {other*xyz}; }
    GEICoord operator*(const double& other)
    { return {other*xyz}; }
    GEICoord operator-()
    { return {-xyz}; }
    // .....................................

    // .....................................
    // Spherical
    double r()      const ; // { return sqrt(xyz.mag2()) ;}
    double alpha()  const ; //  { return RAD_2_DEG*atan2(xyz.y(), xyz.x()) ;}
    double delta()  const ; // { return RAD_2_DEG*asin(xyz.z() / r()) ;}
    // .....................................

    /// Convert to Geocentric Ecliptic coordinates
    virtual GEclipticCoord ToGEcliptic(const Zaki::Physics::Date&) override ;

    /// Convert to Geocentric Ecliptic coordinates
    /// overloaded for UnixTDays input
    virtual GEclipticCoord ToGEcliptic(const double&) override ;
    
    /// Convert to Heliocentric Ecliptic coordinates
    virtual HEclipticCoord ToHEcliptic(const Zaki::Physics::Date&) override ;
    
    /// Convert to Heliocentric Ecliptic coordinates
    /// overloaded for UnixTDays input
    virtual HEclipticCoord ToHEcliptic(const double&) override ;
    
    /// Convert to Geocentric Equatorial coordinates
    virtual GEICoord ToGEI(const Zaki::Physics::Date&) override ;
    
    /// Convert to Geocentric Equatorial coordinates
    /// overloaded for UnixTDays input
    virtual GEICoord ToGEI(const double&) override ;
};

//==============================================================
// The obliquity of the ecliptic (in degrees)
double EclipticObliquity(const Zaki::Physics::Date&) ;
    
// Overloading for UnixTDay input
double EclipticObliquity(const double&) ;


//==============================================================
//                  ostream << overloading
//==============================================================
std::ostream& operator << (std::ostream &, const YPR& );
std::ostream& operator << (std::ostream &, const Coord3D&) ;
std::ostream& operator << (std::ostream &, const GEICoord&) ;
std::ostream& operator << (std::ostream &, const GEclipticCoord&) ;
std::ostream& operator << (std::ostream &, const HEclipticCoord&) ;
//==============================================================

//--------------------------------------------------------------
} // End of namespace Zaki::Physics
//--------------------------------------------------------------

//==============================================================
#endif /*Zaki_Physics_Coordinate_H*/
