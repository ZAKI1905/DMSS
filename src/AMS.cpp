/*
  AMS class

*/

#include <string>

#include "DMSS/AMS.hpp"
//==============================================================
//  AMS Class
// 
// PHYSICAL REVIEW LETTERS 122, 101101 (2019): 
// "The measurement is based on 28.1 × 1e6 electron events
//  collected by AMS from May 19, 2011 to November 12, 2017."  
//      [ 2370 days ]
//--------------------------------------------------------------
// Constructor
AMS::AMS() 
{
  // logger.SetUnit("AMS") ;

  SetName("AMS")              ;
  // SetWrkDir("results/AMS")    ; // This is done by analysis now
  
  NUM_SCALING = 5e-6 ;
  // NUM_SCALING = 8e-7 ;

  
  // Boosted 2017, { 200, 225, 250, 275, 300, 375, 400, 450 }: 7e-4
  // Boosted 2017, { 325, 350, 475, 500, 600, 650 } :  6e-4
  // Boosted 2017, { 425, 550, 700 } : 5e-4
  // Boosted 2017, { 750, 800, 900, 950, 1000, , 1200 } :  4e-4
  // Boosted 2017, { 1150 } : 3e-4 
  // Boosted 2017, { 850, 1050, 1100, 1250, 1300, 1350, 1500, 1600, 1700, 1800, 1900 } : 2e-4
  // Boosted 2017, { 1400, 1900, 2100, 2400, 2600, 2700 } :  9e-5
  // Boosted 2017, { 1700, 2000, 2200, 2500, 3000 } : 8e-5
  // Boosted 2017, { 2300, 2800, 2900 } : 7e-5
  // Boosted 2017, { 3100, 3400, 3500, 3600, 3900 } : 6e-5
  // Boosted 2017, { 3300, 3800, 4000 } : 5e-5
  // Boosted 2017, { 3700, 4200, 4300, 4500, 4900 } : 4e-5
  // Boosted 2017, { 3200, 4400, 4600, 5000 } :  3e-5
  // Boosted 2017, { 4100, 4700, 4800 } :  2e-5


  // Boosted 2028, { 225, 250, 275, 325, 350, 375, 425, 450, 475, 550,
  //                 600, 700, 750, 800, 850, 950, 1500, 4900}:  5e-6
  // Boosted 2028, {200, 300}:  3e-6
  // Boosted 2028, {400}:  1e-6
  // 
  // Boosted 2028, {500, 650, 1000, 1050, 1100, 1150, 1200, 1250, 1300,
  //  1350, 1400, 1600, 1700, 1800, 1900, 2000, 2100,
  //  2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000, 3100, 3200, 
  //  3300, 3400, 3500, 3600, 3700, 3800, 3900, 4000, 4100, 4200, 4300, 
  //  4400, 4500, 4600, 4700, 4800, 5000 } : 8e-7

  SetOrbitInclination(51.6)   ; //  degree
  SetConeFOV(48)              ; // degrees (half-cone angle)
  
  // SetKappa: Overall factor that multiplies the signal numbers.
  // SetKappa(1)                 ;  // should be removed
  
  SetActualTimeDuration({{2011, 5, 19, {0, 0, 0}}, {2017, 11, 13, {0, 0, 0}}}); 
  // SetTimeDuration({{2011, 5, 19, {0, 0, 0}}, {2021, 3, 27, {0, 0, 0}}}) ; 

  SetEccentricity(0)          ;
  SetOrbitPeriod(92.68)       ; // minutes
  SetOrbitHeight(400)         ; // km
  SetInitRightAscNode(-54.6747); //degree
  SetInitTrueAnomaly(0)       ; // theta_0 in degree

  // Table 1 from :
  // "Thermal conditions on the International Space Station: Effects
  //  of operations of the station Main Radiators on the Alpha
  //  Magnetic Spectrometer"
  // Nuclear Instruments and Methods in Physics Research A 815 (2016) 50–56.
  // http://dx.doi.org/10.1016/j.nima.2016.01.048.
  SetYPR(-4, -2, +0.7-12)     ; // Yaw, Pitch, Roll in degree
  FindPrecessionRate()        ;

  Zaki::String::Directory dir(__FILE__) ;
  if(dir.ParentDir().ParentDir().Str() == dir.ParentDir().Str())
    dir = "" ;
  else
    dir = dir.ParentDir().ParentDir() + "/" ;

  ImportData("AMS Electron", dir +"data/AMS_e-.dat", {DMSolarSignal::Electron}) ;

  Z_LOG_NOTE("AMS constructor called for: " + PtrStr()) ;
}

//--------------------------------------------------------------
// Copy constructor
AMS::AMS(const AMS& other) 
{
  *this = other ;
  Z_LOG_NOTE("AMS copy constructor: from " + other.PtrStr() + " --> " + PtrStr()) ;
}

//--------------------------------------------------------------
// Destructor
AMS::~AMS() 
{
  Z_LOG_NOTE("AMS destructor called for: " + PtrStr()) ;
}

//--------------------------------------------------------------
// IClone method
AMS* AMS::IClone() const
{
  Z_LOG_NOTE("AMS clone called for: " + PtrStr()) ;
  return (new AMS(*this)) ;
}

//--------------------------------------------------------------
double AMS::ExpTimeFrac(double energy) const
{
        if ( 0 < energy && energy <= 1 )
    return 0;

  else  if ( 1 < energy && energy <= 2 )
    return 0.06;

  else  if ( 2 < energy && energy <= 35 )
    return (2.242*energy + 1.515)*0.01;

  else  if ( 35 < energy)
    return 0.8;
  else 
    return 0; 
}

//--------------------------------------------------------------
double AMS::Acceptance(double energy) const
{
  double delta = -(0.04 + 0.03)/2 ;
  double A_geom = 550*1e-4 ; // in  units of m^2*sr

        if ( 0 < energy && energy <= 10 )
    return 1 * (1 + delta ) * A_geom;

  else  if ( 10 < energy && energy <= 100 )
    return (0.908 - 0.00078 * energy)* (1 + delta ) * A_geom ;

  else  if ( 100 < energy && energy <= 1400)
    return (0.844 - 0.00014 * energy) * (1 + delta ) * A_geom ;
  else 
    return 0; 
}

//--------------------------------------------------------------

//==============================================================
