/*
  CALET class
  
  Proton:
    "The observation period used in this analysis is from October 13, 2015 to August 31, 2018 (1054 days). "
    https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.122.181102
  
  Electron-Positron:
    " ... electron + positron spectrum from 11 GeV to 4.8 TeV. "
    Using 780 days of flight data from October 13, 2015 to November 30, 2017
    https://arxiv.org/pdf/1806.09728.pdf
    https://link.aps.org/doi/10.1103/PhysRevLett.120.261102

*/

#include "DMSS/CALET.hpp"
//==============================================================
//  CALET Class
//--------------------------------------------------------------
// Constructor
CALET::CALET() 
{
  // logger.SetUnit("CALET")          ;
  
  SetName("CALET")              ;
  // SetWrkDir("results/CALET")    ; // This is done by analysis now
  
  // NUM_SCALING = 1 ;
  NUM_SCALING = 2e-6 ;

  // Boosted (2017): {200} ~ 2e-5
  // Boosted (2017): {225, 250, 300, 350, 375, 400, 425, 450, 475,
                      // 500, 550, 600, 650, 700, 750, 800, 850, 900,
                      // 950, 1000, 1050, 1100, 1150, 1200, 1250,
                      // 1300, 1350, 1400, 1500} ~ 8e-6
                    
  // Boosted (2017): {275, 325, 1600, 1700, 1800, 1900, 2000, 2100, 
                    // 2200, 2300, 2400, 2500, 
                    // 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400,
                    // 3500, 3600, 3700, 3800, 3900, 4000, 4100, 4200,
                    // 4300, 4500, 4600, 4700, 4800, 4900, 5000} ~ 4e-6

  //  Boosted (2017): { 4400 } ~ 3e-6
                    

  // Boosted (2020): {200} ~ 1e-5
  // Boosted (2020): {225, 250, 275, 300, 325, 350, 375, 400, 425,
  //                  450, 475, 500, 550, 600, 650, 700, 
                  //  750, 800, 850, 
                  //  1000, 1050, 1100, 1200, 1250,
                  //  1300, 1350, 1400, 1500, 1600, 1700, 1800, 1900,
                  //  2000, 2100, 2200, 2300, 2400, 
                  //  2600, 2700, 2800, 2900, 3000,
                  //  3100, 3200, 3300, 3500,
                  //  3700, 3800, 3900, 4000  } ~ 4e-6
                  
  // Boosted (2020): { 900, 950, 1150, 2500, 3400, 3600,
                  //   4100, 4200, 4300, 4400, 4500,
                  //   4600, 4700, 4800, 4900 5000} ~ 2e-6


  SetOrbitInclination(51.6)   ; //  degree
  SetConeFOV(45)              ; // degrees (from zenith)
  
  // SetKappa(1)                 ; // should be removed
  
  SetActualTimeDuration({{2015, 10, 13, {0, 0, 0}}, {2017, 11, 31, {0, 0, 0}}}); 
  SetEccentricity(0)          ;
  SetOrbitPeriod(92.68)       ; // minutes
  SetOrbitHeight(400)         ; // km
  SetInitRightAscNode(-168.0292); //degree
  SetInitTrueAnomaly(0)       ; // theta_0 in degree
  SetYPR(-4, -2, +0.7)        ; // Yaw, Pitch, Roll in degree
  FindPrecessionRate()        ;

  Zaki::String::Directory dir(__FILE__) ;
  if(dir.ParentDir().ParentDir().Str() == dir.ParentDir().Str())
    dir = "" ;
  else
    dir = dir.ParentDir().ParentDir() + "/" ;
  ImportData("Calet e+e-", dir + "data/CALET_e+e-.dat", {DMSolarSignal::Electron, DMSolarSignal::Positron}) ;

  Z_LOG_NOTE(("CALET constructor called for: " + PtrStr()).c_str()) ;

}

//--------------------------------------------------------------
// Copy constructor
CALET::CALET(const CALET& other) 
{
  *this = other ;
  Z_LOG_NOTE(("CALET copy constructor: from " + other.PtrStr() + " --> " + PtrStr()).c_str()) ;

}

//--------------------------------------------------------------
// IClone method
CALET* CALET::IClone() const
{
  Z_LOG_NOTE(("CALET::Clone() called for: " + PtrStr()).c_str()) ;

  return (new CALET(*this)) ;
}

//--------------------------------------------------------------
// Destructor
CALET::~CALET() 
{
  Z_LOG_NOTE(("CALET destructor called for: " + PtrStr()).c_str()) ;
}

//--------------------------------------------------------------
// arXiv:1806.09728, Page 3, Bottom_Right:
// "Total live time in this period was 15,811 hours, 
// corresponding to a live time fraction of 84%."
//.......................................
// arXiv:1712.01711, Page 3, Top_Right:
// "On-orbit data collection has been continuous and very stable."
double CALET::ExpTimeFrac(double energy) const
{
  return 0.84 ; 
}

// --------------------------------------------------------------
// arXiv:1806.09728, Page 3, Bottom_Left:
// "... geometrical factor of ∼1040 cm2sr for high-energy electrons."
//.......................................
// arXiv:1712.01711, Page 3, Top_Left:
// " high-energy (>10 GeV).... nearly independent of energy."
//.......................................
// arXiv:1712.01711, Page 4, Top_Left:
// "Combined efficiency of preselection for electrons is very high:
//  >90% above 30 GeV to 3 TeV 
//  85%  at    20 GeV at variance with only
//  60%  at    10 GeV due to lower trigger efficiency."
//.......................................
// arXiv:1712.01711, Page 6, Top_Left:
// ε(E) is the detection efficiency for electrons defined as 
// the product of trigger, preselection, track reconstruction
//  and electron identification efficiencies..."
//
// "Based on the MC simulations, the total efficiency is very 
// stable with energy up to 3 TeV: 73% ± 2%.""
double CALET::Acceptance(double energy) const
{
  double A_geom = 1040*1e-4 ; // in  units of m^2*sr

  double eff = 0.73 ; // efficiency 
  
  return A_geom*eff; 
}

//--------------------------------------------------------------

//==============================================================
