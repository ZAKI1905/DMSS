/*
  DAMPE class

*/

#include "DMSS/DAMPE.hpp"

//==============================================================
//  DAMPE Class
// https://arxiv.org/pdf/1711.10981.pdf, Page 3, 2nd paragraph
// "The results reported in this work are based on data recorded
// between December 27, 2015 and June 8, 2017." [ ~530 days ]
// 
//--------------------------------------------------------------
// Constructor
DAMPE::DAMPE() 
{
  // logger.SetUnit("DAMPE") ;
  SetName("DAMPE") ;
  // SetWrkDir("results/DAMPE")      ; // This is done by analysis now

  NUM_SCALING = 1  ;
  // NUM_SCALING = 1e-6 ;

  // Boosted (2017): {200, 225, 250, 275, 300, 325, 350, 375,
                  //  425, 450, 475, 500, 550, 600, 650} ~ 2e-5
  // Boosted (2017): {400} ~ 5e-5
  // Boosted (2017): { 700, 750, 800, 850 } ~ 1e-6
  // Boosted (2017): {900, 950, 1000, 1500, 1600, 1800, 1900} ~ 5e-7
  // Boosted (2017): {1050, 1100, 1150, 1200, 1250, 1300, 1350, 1400} ~ 1e-7
  // Boosted (2017): {1700, 2000, 2100, 2200, 2400, 2500, 2600,
                  //  2700, 2800, 2900, 3400, 3500, 3600, 3700,
                  //  3900, 4400, 4500, 4600, 4800} ~ 5e-8
  // Boosted (2017): { 3000, 4900 } ~ 4e-8 
  // Boosted (2017): { 2300, 3100, 3200, 3300, 3800, 4000, 4100,
                  //  4200, 4300 } ~ 1e-8 
  // Boosted (2017): { 4700 } ~ 8e-8
  // Boosted (2017): { 5000 } ~ 6e-8



  // Boosted (2020): {200, 225, 250, 275, 300, 325, 350, 375,
                  //  400, 425, 475, 500} ~ 2e-5
  // Boosted (2020): {450} ~ 5e-5
  // Boosted (2020): {550, 600, 650, 700, 750, 800, 850, 900,
                  //  950, 1000, 1050, 1100, 1150, 1200, 1250,
                  // 1300, 1350, 1400, 1600, 1700, 1800, 1900, 2000 } ~ 1e-5
  // Boosted (2020): {2400, 3000, 3200} ~ 9e-6
  // Boosted (2020): {1500, 2100, 2200, 2300, 2500, 2600, 2700,
                  //  2800, 3100, 3400, 3600, 3700, 3800, 3900, 4000} ~ 5e-6
  // Boosted (2020): {2900, 3500, 4200} ~ 1e-6
  // Boosted (2020): { 3300 } ~ ???
  // Boosted (2020): {4100, 4400, 4500, 4600} ~ 8e-6
  // Boosted (2020): {4300, 4700, 4800, 4900, 5000} ~ 4e-6


  // SetConeFOV(62)                  ; // degrees (half-cone angle)

  // SetKappa(1)                     ; // should be removed
  
  SetOrbitInclination(97.406)     ; //  degree
  SetActualTimeDuration({{2015, 12, 27, {0, 0, 0}}, {2017, 6, 9, {0, 0, 0}}}); 
  SetEccentricity(0)              ;
  SetOrbitPeriod(95)              ; // minutes
  SetOrbitHeight(500)             ; // km
  SetInitRightAscNode(5.326)      ; //degree
  SetInitTrueAnomaly(0)           ; // theta_0 in degree
  SetYPR(0, 0, 0)                 ; // Yaw, Pitch, Roll in degree
  
  // Forcing the orbit to precess exactly with Sun
  ForceSunSync()                  ;
  // FindPrecessionRate()            ;

  //  DAMPE energy bins taken from:
  //  https://www.nature.com/articles/nature24475
  Zaki::String::Directory dir(__FILE__) ;
  if(dir.ParentDir().ParentDir().Str() == dir.ParentDir().Str())
    dir = "" ;
  else
    dir = dir.ParentDir().ParentDir() + "/" ;
  ImportData("DAMPE e+e-", dir +"data/DAMPE_e+e-.dat",  {DMSolarSignal::Electron, DMSolarSignal::Positron});

  Z_LOG_NOTE(("DAMPE constructor called for: " + PtrStr()).c_str()) ;

}

//--------------------------------------------------------------
// Copy constructor
DAMPE::DAMPE(const DAMPE& other) 
{
  *this = other ;
  Z_LOG_NOTE(("DAMPE copy constructor: from " + other.PtrStr() + " --> " + PtrStr()).c_str()) ;

}

//--------------------------------------------------------------
// Clone method
DAMPE* DAMPE::IClone() const
{
  Z_LOG_NOTE(("DAMPE::Clone() called for: " + PtrStr()).c_str()) ;

  return (new DAMPE(*this)) ;
}

//--------------------------------------------------------------
// Destructor
DAMPE::~DAMPE() 
{
  Z_LOG_NOTE(("DAMPE destructor called for: " + PtrStr()).c_str()) ;
}

//--------------------------------------------------------------
// Acceptance in units of m^2*sr, taken from:
// Table 1 in https://www.nature.com/articles/nature24475
double DAMPE::Acceptance(double energy) const
{
        if ( 24.0 < energy && energy <= 27.5 )
    return 0.256 ;

  else  if ( 27.5 < energy && energy <= 31.6 )
    return 0.259 ;

  else  if ( 31.6 < energy && energy <= 36.3)
    return 0.261 ;

  else  if ( 36.3 < energy && energy <= 41.7)
    return 0.264 ;

  else  if ( 41.7 < energy && energy <= 47.9)
    return 0.266 ;

  else  if ( 47.9 < energy && energy <= 55.0)
    return 0.269 ;

  else  if ( 55.0 < energy && energy <= 63.1)
    return 0.272 ;
  
  else  if ( 63.1 < energy && energy <= 72.4)
    return 0.275 ;
  
  else  if ( 72.4 < energy && energy <= 83.2)
    return 0.277 ;
  
  else  if ( 83.2 < energy && energy <= 95.5)
    return 0.279 ;
  
  else  if ( 95.5 < energy && energy <= 109.7)
    return 0.283 ;
  
  else  if ( 109.7 < energy && energy <= 125.9)
    return 0.282 ;
  
  else  if ( 125.9 < energy && energy <= 144.5)
    return 0.286 ;
  
  else  if ( 144.5 < energy && energy <= 166.0)
    return 0.287 ;
  
  else  if ( 166.0 < energy && energy <= 190.6)
    return 0.288 ;

  else  if ( 190.6 < energy && energy <= 218.8)
    return 0.288 ;

  else  if ( 218.8 < energy && energy <= 251.2)
    return 0.290 ;
  
  else  if ( 251.2 < energy && energy <= 288.4)
    return 0.291 ;
  
  else  if ( 288.4 < energy && energy <= 331.1)
    return 0.291 ;
  
  else  if ( 331.1 < energy && energy <= 380.2)
    return 0.290 ;
  
  else  if ( 380.2 < energy && energy <= 436.5)
    return 0.292 ;
  
  else  if ( 436.5 < energy && energy <= 501.2)
    return 0.291 ;
  
  else  if ( 501.2 < energy && energy <= 575.4)
    return 0.289 ;
  
  else  if ( 575.4 < energy && energy <= 660.7)
    return 0.288 ;
  
  else  if ( 660.7 < energy && energy <= 758.6)
    return 0.285 ;
 
  else  if ( 758.6 < energy && energy <= 871.0)
    return 0.284 ;
  
  else  if ( 871.0 < energy && energy <= 1000.0)
    return 0.278 ;
  
  else  if (1000.0 < energy && energy <= 1148.2)
    return 0.276 ;

  else  if (1148.2 < energy && energy <= 1318.3)
    return 0.274 ;

  else  if (1318.3 < energy && energy <= 1513.6)
    return 0.267 ;
  
  else  if (1513.6 < energy && energy <= 1737.8)
    return 0.263 ;
  
  else  if (1737.8 < energy && energy <= 1995.3)
    return 0.255 ;
  
  else  if (1995.3 < energy && energy <= 2290.9)
    return 0.249 ;
  
  else  if (2290.9 < energy && energy <= 2630.3)
    return 0.243 ;
  
  else  if (2630.3 < energy && energy <= 3019.9)
    return 0.233 ;
  
  else  if (3019.9 < energy && energy <= 3467.4)
    return 0.227 ;
  
  else  if (3467.4 < energy && energy <= 3981.1)
    return 0.218 ;
  
  else  if (3981.1 < energy && energy <= 4570.9)
    return 0.210 ;

  else 
    return 0; 
}

//--------------------------------------------------------------
//==============================================================
