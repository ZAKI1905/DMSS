#ifndef Zaki_Physics_Constants_H
#define Zaki_Physics_Constants_H

#include <math.h>
#include <vector>
#include <algorithm>
#include <string>

//--------------------------------------------------------------
namespace Zaki::Physics
{

//......................................
// Forward declrations
//......................................
struct Date;
//==============================================================
//                      Unit Conversion
//==============================================================
//......................................
//          Length
//......................................
extern const double  M_2_CM     ;           // meter to cm conversion
extern const double  KM_2_M     ;           // kilometer to meter
extern const double  M_2_KM     ;           // meter to kilometer
extern const double  CM_2_M     ;           // cm to meter conversion
extern const double  CM_2_GEV   ;           // cm to GeV conversion
extern const double  GEV_2_CM   ;           // GeV to cm conversion
extern const double  M_2_GEV    ;           // m to GeV conversion
extern const double  KM_2_GEV   ;           // kilometer to GeV
extern const double  GEV_2_M    ;           // GeV to m conversion
extern const double  GEV_2_KM   ;           // GeV to km
extern const double  AU_2_KM    ;           // AU to km
extern const double  AU_2_GEV   ;           // AU to GeV

//......................................
//          Time
//......................................
extern const double  SEC_2_GEV  ;           // sec to GeV^-1 conversion
extern const double  GEV_2_S    ;           // GeV to sec conversion
extern const double  GYR_2_YR   ;           // Giga year to year
extern const double  YR_2_DAY   ;           // Year to days
extern const double  DAY_2_SEC  ;           // Year to days
extern const double  YR_2_SEC   ;           // Year to days
extern const double  GYR_2_SEC  ;           // Giga year to seconds
extern const double  GYR_2_GEV  ;           // Giga year to GeV^-1
extern const double  MIN_2_DAY  ;           // minutes to day

//......................................
//          Temperature
//......................................
extern const double  KEL_2_GEV  ;           // Kelvin to GeV conversion
extern const double  GEV_2_KEL  ;           // GeV to Kelvin conversion

//......................................
//           Mass
//......................................
extern const double  KG_2_GR    ;           // Kilogram to gram conversion
extern const double  GR_2_KG    ;           // gram to Kilogram conversion
extern const double  GEV_2_GR   ;           // GeV to gram conversion
extern const double  GR_2_GEV   ;           // gram to GeV conversion
extern const double  KG_2_GEV   ;           // kilogram to GeV conversion
extern const double  GEV_2_KG   ;           // GeV to kilogram conversion
extern const double  AMU_2_GEV  ;           // amu to GeV conversion
extern const double  GEV_2_AMU  ;           // GeV to amu conversion

//......................................
//           Angle
//......................................
extern const double RAD_2_DEG    ;          // Radian to degree
extern const double DEG_2_RAD    ;          // Degree to radian

//......................................
//    Compact Star Unit Conversions
//......................................
extern const double MEV_2_FM        ;       // MeV to 1/fm
extern const double MEV_FM3_2_G_CM3 ;       // MeV/fm^3 to g/cm^3
extern const double G_CM3_2_Dyn_CM2 ;       // g/cm^3 to dyne/cm^2

extern const double MEV_FM3_2_Dyn_CM2 ;
extern const double FM4_2_G_CM3       ;     // 1/fm^4 to g/cm^3
extern const double FM4_2_Dyn_CM2     ;     // 1/fm^4 to dyne/cm^2

//==============================================================
//                      Constants
//==============================================================
//......................................
//              Physical constants
//......................................
extern const double  ALPHA_EM       ;       // EM structure const doubleant
extern const double  Q_E            ;       // Electric charge
extern const double  Q_E_SQRD_MeVfm ;       // Electric charge in MeV.fm

extern const double NEWTON_G_SI     ;       //  m^3/(kg s^2)
extern const double NEWTON_G_GEV    ;       //  in 1/GeV^2

extern const double LIGHT_C_M_S     ; // Speed of light in m/s
extern const double LIGHT_C_KM_S    ; // Speed of light in km/s

//......................................
//              Particles
//......................................
// Electron
extern const double ELECTRON_M_KG   ;       // Mass of the electron in kilograms
extern const double ELECTRON_M_GEV  ;       // Mass of the electron in GeV
extern const double ELECTRON_M_FM   ;       // Mass of the electron in 1/fm

// Muon
extern const double MUON_M_GEV      ;       // GeV
extern const double MUON_M_MEV      ;       // MeV
extern const double MUON_M_FM       ;       // fm^-1

// Tau Lepton
extern const double TAU_M_GEV      ;       // GeV
extern const double TAU_M_MEV      ;       // MeV
extern const double TAU_M_FM       ;       // fm^-1

// Proton
extern const double PROTON_M_GEV    ;       // GeV
extern const double PROTON_M_MEV    ;       // MeV
extern const double PROTON_M_FM     ;       // fm^-1

// Neutron
extern const double NEUTRON_M_GEV   ;       // GeV
extern const double NEUTRON_M_MEV   ;       // MeV
extern const double NEUTRON_M_FM    ;       // fm^-1

// Sigma^- Baryon
extern const double SIGMA_MINUS_M_GEV   ;       // GeV
extern const double SIGMA_MINUS_M_MEV   ;       // MeV
extern const double SIGMA_MINUS_M_FM    ;       // fm^-1

// Sigma^0 Baryon
extern const double SIGMA_ZERO_M_GEV   ;       // GeV
extern const double SIGMA_ZERO_M_MEV   ;       // MeV
extern const double SIGMA_ZERO_M_FM    ;       // fm^-1

// Sigma^+ Baryon
extern const double SIGMA_PLUS_M_GEV   ;       // GeV
extern const double SIGMA_PLUS_M_MEV   ;       // MeV
extern const double SIGMA_PLUS_M_FM    ;       // fm^-1

// Lambda^0 Baryon
extern const double LAMBDA_ZERO_M_GEV   ;       // GeV
extern const double LAMBDA_ZERO_M_MEV   ;       // MeV
extern const double LAMBDA_ZERO_M_FM    ;       // fm^-1

// Xsi^- Baryon
extern const double XSI_MINUS_M_GEV   ;       // GeV
extern const double XSI_MINUS_M_MEV   ;       // MeV
extern const double XSI_MINUS_M_FM    ;       // fm^-1

// Xsi^0 Baryon
extern const double XSI_ZERO_M_GEV   ;       // GeV
extern const double XSI_ZERO_M_MEV   ;       // MeV
extern const double XSI_ZERO_M_FM    ;       // fm^-1

//......................................
//              Sun
//......................................
// Sun radius 
extern const double SUN_R_CM        ;       // in cm
extern const double SUN_R_KM        ;       // in km
extern const double SUN_R_GEV       ;       // in GeV

// Sun core temperature
extern const double  SUN_T_KEL      ;       // Kelvin
extern const double  SUN_T_GEV      ;       //  GeV

// Sun core density
extern const double  SUN_RHO_SI     ;       // g / cm^3
extern const double  SUN_RHO_GEV    ;       //  GeV^4

// Sun's age
extern const double  SUN_AGE_GYR    ;       // Sun's age in Gyr
extern const double  SUN_AGE_GEV    ;       // Sun's age in GeV

// Solar precession rate
extern const double SOLAR_PREC_RATE ;       // degree / day  

//......................................
//              Dark Matter
//......................................
// Dark matter density in the Solar system
extern const double  DM_RHO       ;         // in GeV/cm^3
extern const double  DM_RHO_GEV   ;         // in GeV^4

// cross section needed for relic abundance (cm^3/s)
extern const double SIG_V_REL     ;

// Angular average of the DM velocity distribution with respect to the solar rest
// frame neglecting the gravitational attraction of the Sun.
extern const double  F_Sun        ;         // in natural units 

//......................................
//              Earth
//......................................
extern const double EARTH_ROT_AXIS  ;       // Degrees
extern const double EARTH_RADIUS    ;       // km

// Earth' s second dynamic form factor
extern const double EARTH_J2        ;       // Earth's J2 factor

// Distance betwen the Sun & Earth
extern const double EARTH_2_SUN     ;       // km
extern const double EARTH_2_SUN_GEV ;       // GeV

//......................................
//              Date & Time
//......................................
// Noon on January 1, 2000
extern const Date J2000 ;   

// 00:00:00 UTC on 1 January 1970
extern const Date UnixEpoch ;   

// 00:00:00 UTC on 1 January 1995 (for root)
extern const Date JAN1st1995 ; 

//...........................................................
//            Vernal equinox (UT)
// Ref:
//  1. https://en.wikipedia.org/wiki/March_equinox
//  2. http://www.astropixels.com/ephemeris/soleq2001.html
//...........................................................
// Vernal equinox in 2001 (in UT)
extern const Date VERNAL_2001;

// Vernal equinox in 2002 (in UT)
extern const Date VERNAL_2002;

// Vernal equinox in 2003 (in UT)
extern const Date VERNAL_2003;

// Vernal equinox in 2004 (in UT)
extern const Date VERNAL_2004;

// Vernal equinox in 2005 (in UT)
extern const Date VERNAL_2005;

// Vernal equinox in 2006 (in UT)
extern const Date VERNAL_2006;

// Vernal equinox in 2007 (in UT)
extern const Date VERNAL_2007;

// Vernal equinox in 2008 (in UT)
extern const Date VERNAL_2008;

// Vernal equinox in 2009 (in UT)
extern const Date VERNAL_2009;

// Vernal equinox in 2010 (in UT)
extern const Date VERNAL_2010;

// Vernal equinox in 2011 (in UT)
extern const Date VERNAL_2011;

// Vernal equinox in 2012 (in UT)
extern const Date VERNAL_2012;

// Vernal equinox in 2013 (in UT)
extern const Date VERNAL_2013;

// Vernal equinox in 2014 (in UT)
extern const Date VERNAL_2014;

// Vernal equinox in 2015 (in UT)
extern const Date VERNAL_2015;

// Vernal equinox in 2016 (in UT)
extern const Date VERNAL_2016;

// Vernal equinox in 2017 (in UT)
extern const Date VERNAL_2017;

// Vernal equinox in 2018 (in UT)
extern const Date VERNAL_2018;

// Vernal equinox in 2019 (in UT)
extern const Date VERNAL_2019;

// Vernal equinox in 2020 (in UT)
extern const Date VERNAL_2020;
//...............................................

//==============================================================

//==============================================================
// Chemical Elements
struct Element
{
    std::string sym = "" ;
    double Z = 0 ; // Atomic number
    double A = 0 ; // Mass number, Standard atomic weight 
    double I = 0 ; // DM capture rate in the Sun
    double v2 = 0; // Typical escape velocity squared in the Sun

    Element() {} 
    Element(const std::string& in_sym, double in_Z, double in_A)
        : sym(in_sym), Z(in_Z), A(in_A) {}
    Element(const std::string& in_sym, double in_Z, double in_A, double in_I, double in_v2)
        : sym(in_sym), Z(in_Z), A(in_A), I(in_I), v2(in_v2) {}

    double mu() { return A*AMU_2_GEV ;} // Molar mass in GeV
};
//==============================================================

//--------------------------------------------------------------
} // End of namespace Zaki::Physics
//--------------------------------------------------------------

#endif /*Zaki_Physics_Constants_H*/
