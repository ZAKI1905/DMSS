/*
  Satellite class

*/

#include <gsl/gsl_integration.h>

// Creating directory
#include <sys/stat.h>
#include <array>

// Root
#include <TDatime.h>
#include <TStopwatch.h>
#include <TF1.h>
#include <TF3.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLine.h>
#include <TGraph2D.h>
#include <TMultiGraph.h>

#include <Zaki/File/CSVIterator.hpp>
#include <Zaki/Vector/Vector_Basic.hpp>
#include <Zaki/String/TextBox.hpp> // For Print()
#include <Zaki/Physics/Coordinate.hpp>
#include <Zaki/Math/GSLFuncWrapper.hpp>

// Local headers
#include "DMSS/Satellite.hpp"

using namespace Zaki::Physics ;
//==============================================================

//--------------------------------------------------------------
// Constructor
Satellite::Satellite() 
  : Prog("Satellite", true)
{}

//--------------------------------------------------------------
Satellite::~Satellite() {}

//--------------------------------------------------------------
/*  

Position of the Sun (in km) as a function of time

Taken from:
https://en.wikipedia.org/wiki/Position_of_the_Sun :

 " These equations, from the Astronomical Almanac,[1][2] can be
  used to calculate the apparent coordinates of the Sun, mean
  equinox and ecliptic of date, to a precision of about 0.01 (36"),
  for dates between 1950 and 2050."

Ref:
1. U.S. Naval Observatory; U.K. Hydrographic Office, 
 H.M. Nautical Almanac Office (2008). 
 "The Astronomical Almanac for the Year 2010". U.S. Govt. 
  Printing Office. p. C5. ISBN 978-0-7077-4082-9.

2. Much the same set of equations, covering the years 1800 to 2200,
  can be found at Approximate Solar Coordinates:
  http://aa.usno.navy.mil/faq/docs/SunApprox.php
  , at the U.S. Naval Observatory website Archived 2016-01-31 at the Wayback Machine.
  https://web.archive.org/web/20160131231447/http://www.usno.navy.mil/USNO
  Graphs of the error of these equations, compared to an accurate
  ephemeris, can also be viewed.

3. ``CALCULATION OF THE SUN, MOON AND ISS POSITIONS''
    NICOLA TOMASSETTI, INFN - Sezione di Perugia, 06122 Perugia, Italy.

  Input should be in days from Unix epoch (1970)
  Output unit is in km
*/
Zaki::Physics::GEICoord Satellite::GetSunPos(double in_date) const
{
  // // Julian day number at noon on January 1, 2000:
  // double J2000 = 2451545 ;

  double days_since_J2000 =  in_date - J2000.UnixTDay() ; 

  // Julian centuries from the epoch J2000
  double t_J2000 =  days_since_J2000 / 36525 ;

  Zaki::Physics::Date tmp_1950 = {1950, 1, 1};
  Zaki::Physics::Date tmp_2050 = {2050, 1, 1};
  if (in_date < tmp_1950.UnixTDay() || in_date > tmp_2050.UnixTDay())
  {
    Z_LOG_ERROR("Year must be in [1950, 2050] range, otherwise the answer is not accurate!") ;
  } 

  // mean longitude (deg)
  // double VL = fmod(280.4665 + 0.9856474*days_since_J2000, 360.0) ;
  double VL = fmod(280.4665 + 36000.76983*t_J2000 + 0.0003032*pow(t_J2000,2),
                   360.0) ;

  //  mean anomaly (rad)
  // double g = fmod(357.52911 + 0.9856003*days_since_J2000, 360.0) * DEG_2_RAD ;
  double M = fmod(357.52911 + 35999.05029*t_J2000 - 0.0001537*pow(t_J2000,2),
                  360.0) * DEG_2_RAD ;

  // The eccentricity of the Earth’s orbit:
  double ecc = 0.016708634 - 0.000042037*t_J2000 - 0.0000001267*pow(t_J2000,2) ;

  // Sun’s equation of the center C (degrees):
  double C = (1.914602 - 0.004817*t_J2000 - 0.000014*pow(t_J2000,2) )*sin(M)
            + (0.019993 - 0.000101*t_J2000 )*sin(2*M)
            + 0.000289*sin(3*M) ;

  // Sun’s geometric longitude (deg):
  double geom_long = VL + C;

  // Sun’s true anomaly (deg):
  double true_anom = M + C;

  // Apparent LONGITUDE ALONG ECLIPTIC (rad)
  // double lambda = (VL + 1.915*sin(g) + 0.020 * sin(2*g)) * DEG_2_RAD;
  
  // Degrees
  double omega = 125.04452 - 1934.136261*t_J2000;
  double lambda = geom_long - 0.00569 - 0.00478*sin(omega*DEG_2_RAD);
  lambda       *= DEG_2_RAD ;

  // obliquity of the ecliptic (rad)
  // double epsilon = (23.439 - 0.0000004*days_since_J2000) * DEG_2_RAD ;
  double epsilon = Zaki::Physics::EclipticObliquity(in_date)*DEG_2_RAD ;

  // APPARENT DECLINATION (rad)
  double delta = asin(sin(epsilon)*sin(lambda)) ;

  // APPARENT RIGHT ASCENSION (rad)
  double alpha = atan2(cos(epsilon)*sin(lambda), cos(lambda)) ;
  
  // std::cout << "\n DECLINATION: " << RAD_2_DEG*delta <<
  //             ", RIGHT ASCENSION: " << RAD_2_DEG*alpha << "\n" ;

  // Distance of the Sun from the Earth, in astronomical units
  // double r_sun = 1.00014 - 0.01671*cos(g) - 0.00014*cos(2*g) ;
  // in AU
  double r_sun = 1.000001018*(1-pow(ecc, 2)) / (1+ecc*cos(true_anom*DEG_2_RAD)) ;

  r_sun *= AU_2_KM ; // in km

  Zaki::Physics::GEICoord xyz = {{ 
    r_sun*cos(alpha)*cos(delta),
    r_sun*sin(alpha)*cos(delta), 
    r_sun*sin(delta)}} ;

  return xyz ;
}

//--------------------------------------------------------------
//  The ecliptic longitude of the Sun (in degrees)
double Satellite::GetSunEclipticLong(double in_Jdate) const
{
  double days_since_J2000 =  in_Jdate - J2000.UnixTDay() ; 

  Zaki::Physics::Date tmp_1950 = {1950, 1, 1};
  Zaki::Physics::Date tmp_2050 = {2050, 1, 1};
  if (in_Jdate < tmp_1950.UnixTDay() || in_Jdate > tmp_2050.UnixTDay())
  {
    Z_LOG_ERROR("Year must be in [1950, 2050] range, otherwise the answer is not accurate!") ;
  } 

  // mean longitude (deg)
  double VL = fmod(280.460 + 0.9856474*days_since_J2000, 360.0) ;

  //  mean anomaly (rad)
  double g = fmod(357.528 + 0.9856003*days_since_J2000, 360.0) * DEG_2_RAD ;

  // LONGITUDE ALONG ECLIPTIC (deg)
  double lambda = (VL + 1.915*sin(g) + 0.020 * sin(2*g));

  return lambda - 360*floor( lambda / 360 ) ;
}

//--------------------------------------------------------------
Zaki::Physics::GEICoord Satellite::GetSatPos(double in_Jdate) const
{
  Zaki::Physics::GEICoord tmp_pos = {{1, 0, 0}} ;

  if (in_Jdate < GetTimeDuration().start.UnixTDay() ||
      GetTimeDuration().end.UnixTDay() < in_Jdate)
  {
    Z_LOG_WARNING("The date is not within the satellite time range.") ;
  }

  // double tmp_time_dur = in_Jdate - GetTimeDuration().start.UnixTDay() ;

  tmp_pos = GetOrbitRadius()*(ROOT::Math::RotationZ(GetRAAN(in_Jdate)*DEG_2_RAD)*
                ROOT::Math::RotationX(GetOrbitInclination()*DEG_2_RAD)*
                ROOT::Math::RotationZ(GetTrueAnomaly(in_Jdate)*DEG_2_RAD)*
                tmp_pos.XYZ()) ;
  
  return tmp_pos ;
}

//--------------------------------------------------------------
// Tangent to the orbit at satellite position in GEI coordinates
Zaki::Physics::GEICoord Satellite::GetSatTan(double in_Jdate) const
{
  Zaki::Physics::GEICoord tmp_pos = {{0, 1, 0}} ;

  if (in_Jdate < GetTimeDuration().start.UnixTDay() ||
      GetTimeDuration().end.UnixTDay() < in_Jdate)
  {
    Z_LOG_WARNING("The date is not within the satellite time range.") ;
  }

  // double tmp_time_dur = in_Jdate - GetTimeDuration().start.UnixTDay() ;

  tmp_pos = ROOT::Math::RotationZ(GetRAAN(in_Jdate)*DEG_2_RAD)*
                ROOT::Math::RotationX(GetOrbitInclination()*DEG_2_RAD)*
                ROOT::Math::RotationZ(GetTrueAnomaly(in_Jdate)*DEG_2_RAD)*
                tmp_pos.XYZ() ;
  
  return tmp_pos ;
}

//--------------------------------------------------------------
// Overriding the base method
void Satellite::SetWrkDir(const Zaki::String::Directory& input) 
{
  wrk_dir = input + "/" + name ;
  set_wrk_dir_flag = true ;

  // ............ Creating a directory ............
  if (mkdir(wrk_dir.Str().c_str(), ACCESSPERMS) == -1) 
  {
    Z_LOG_NOTE("Directory '"+input.Str()+"' wasn't created, because: "+strerror(errno)+".") ; 
  }
  else
    Z_LOG_INFO(("Directory '" + wrk_dir.Str() + "' created.").c_str()); 
  // .................................................

  results.SetWrkDir(wrk_dir) ;
}

//--------------------------------------------------------------
std::string Satellite::GetOrbitType() const 
{
  double i = GetOrbitInclination() ; // degree

  if( 90 < i && i < 180)
    return "Retrograde" ;
  else if ( 0 < i && i < 90 )
    return "Prograde" ;
  else if ( i == 0 || i == 180)
    return "Equatorial" ;
  else // i == 90
    return "Polar";
  
}

//--------------------------------------------------------------
// Set Filed-Of-View (FOV) defined as the half-angle from the normal
// to the surface of the detector and in a conical shape.
void Satellite::SetConeFOV(const double& fov_in)             
{
  ConeFOV = fov_in ;
  set_cone_fov_flag = true ;
}

//--------------------------------------------------------------
// Get Filed-Of-View (FOV) defined as the half-angle from the normal
// to the surface of the detector and in a conical shape.
double Satellite::GetConeFOV(std::string in_unit) const            
{
  if( in_unit == "deg")
    return ConeFOV ;
  else if ( in_unit == "sr")
    return 2*M_PI*(1-cos(ConeFOV*DEG_2_RAD)) ;
  else
  {
    Z_LOG_ERROR("Unit must be in 'sr' or 'deg'   -->   assuming 'deg'.") ;
    return ConeFOV ;
  }
}

//--------------------------------------------------------------
void Satellite::Print() const
{
  using namespace Zaki::String;

  TextBox text_box
  (
    {
      "Satellite Name: " + GetName(), "\n",
      "Orbit Inclination: "+ ToString(GetOrbitInclination()) + "\u00B0", "\n",
      "Orbit Type: " + ToString(GetOrbitType()), "\n",
      "Eccentricity: " + ToString(GetEccentricity()), "\n",
      "Time Duration: " + ToString(GetTimeDuration()), "\n",
      "Initial Right Ascention of the Ascending Node (RAAN):",
      "Angle:  " + ToString(GetInitRightAscNode()) + "\u00B0,    Local Time:  " 
        + ToString(GetLocalTAscNode(GetTimeDuration().start.UnixTDay())), "\n",
      "Orbit Precession Rate: " + ToString(GetPrecessionRate()) + "\u00B0/day", "\n",
      "Orbit Height: " + ToString(GetOrbitHeight()) + " km,   Orbit Radius: " 
        + ToString(GetOrbitRadius()) + " km", "\n",
      "Orbit Period: " + ToString(GetOrbitPeriod()) + " min,   Orbit Frequency: " 
        + ToString(GetOrbitAngFreq()) + " rad/min", "\n",
      "Orientation:",
      ToString(GetYPR()), "\n",
      "Conical FOV: " + ToString(GetConeFOV()) + "\u00B0"
    }
  ) ;
  text_box.SetAlignment(TextBox::center) ;
  text_box.Print() ;
}

//--------------------------------------------------------------
// Checks if the satelite is in Earth's shadow
// Ref: "Flight and Orbital Mechanics"
// https://ocw.tudelft.nl/courses/flight-orbital-mechanics/
//  "Methods of Orbit Determination" [Escobal, 1976]
bool Satellite::IsInShadow(const double& t_JD) 
{
  // cosine of the angle between the vector from the Sun to 
  // Earth, and the vector from Earth to the satelite
  double cs = -GetSunPos(t_JD).XYZ().Dot(GetSatPos(t_JD).XYZ()) ;

  // * 1. Satellite on night-side of the Earth
  if (cs <= 0 )  return false ;


  // Normalizing so that cs = cosine
  cs       *= 1.0 / sqrt(GetSunPos(t_JD).XYZ().mag2()) ;
  cs       *= 1.0 / sqrt(GetSatPos(t_JD).XYZ().mag2()) ;


 
  double a = sqrt(GetSatPos(t_JD).XYZ().mag2() )*sin(acos(cs));

  // * 2. Satellite hides behind the Earth
  bool cond_2 = a < EARTH_RADIUS ;

  return  cond_2 ;
}

//--------------------------------------------------------------
void Satellite::ImportExposure(const Zaki::String::Directory& f_name)
{
  std::ifstream     file( (wrk_dir + "/" + f_name).Str());

  // Error opening the file
  if (file.fail()) 
  {
    Z_LOG_ERROR("File '"+(wrk_dir + "/" + f_name).Str() +"' cannot be opened!") ;
    Z_LOG_ERROR("Importing exposure failed!") ;
    exit(EXIT_FAILURE) ;
    return ;
  }

  // Reading the input file
  for(Zaki::File::CSVIterator loop(file); loop != Zaki::File::CSVIterator(); ++loop)
  {
    exposure_set.push_back(std::stof((*loop)[0])) ;
  }

  Zaki::Math::Range<size_t> t_range  = {0, exposure_set.size()};
  unsigned int bin_num   = exposure_set.size() ;

  TH1F tmp_exp_hist(("exp_"+GetName()).c_str(), 
                    ("Exposure_"+GetName()).c_str(), 
                    bin_num, t_range.min, t_range.max);

  // tmp_exp_hist.SetBinContent(0, 0) ; // Underflow
  // Filling the Exposure histogram
  for(size_t i=0; i<exposure_set.size(); i++)
  {
    tmp_exp_hist.SetBinContent(i, exposure_set[i])  ;
  }

  exp_hist = tmp_exp_hist;

  Z_LOG_INFO("Exposure data imported from: "+ (wrk_dir+f_name).Str()+".") ;
  set_exposure_eval_flag = true ;
}

//--------------------------------------------------------------
std::vector<double> Satellite::GetExposure() const
{
  return exposure_set;
} 

//--------------------------------------------------------------
void Satellite::PlotExposure(const Zaki::String::Directory& f_name, double energy)
{

  PlotExposureTime(f_name, energy) ;

  Z_LOG_INFO("Plotting the total exposure vs. energy...") ;
  PlotExposureEnergy(f_name) ;
}

//--------------------------------------------------------------
//  Making plots using root
// Wrapper for: Satellite::PlotExposureEnergy
double Satellite::ExpEnergyFuncWrapper(double *x, double *par)
{

  double exp_tot = exp_hist.Integral() ;
  exp_tot       *= 24*3600*ExpTimeFrac(*x) ;
  exp_tot       *= Acceptance(*x) / GetConeFOV("sr");

  return exp_tot;
}

//--------------------------------------------------------------
void Satellite::PlotExposureEnergy(const Zaki::String::Directory& f_name)
{
  if(!set_exposure_eval_flag)
  {
    Z_LOG_ERROR("Plot failed because exposure hasn't been evaluated!");
    exit(EXIT_FAILURE) ;
    return ;
  }

  // if(!set_energy_range_flag)
  // {
  //   Z_LOG_ERROR("Plot failed because energy range hasn't been set!");
  //   return ;
  // }

  Z_LOG_INFO("Plotting the total exposure vs. energy...") ;

  double e_min = GetData()[0].GetERange().min ;
  double e_max = GetData()[GetData().size()-1].GetERange().max ;

  TF1 ExpPlot("Exposure(e)", this, &Satellite::ExpEnergyFuncWrapper,
              e_min, e_max, 0);

  std::string date_range_str ;
  date_range_str = GetTimeDuration().start.StrForm("#Y-#M-#D")
                  + "   to   "
                  + GetTimeDuration().end.StrForm("#Y-#M-#D");

  char tmp[200] ;
  sprintf(tmp, "#splitline{Total Exposure vs. Energy}{%s}; E (GeV); Exposure (m^2.s)",
          date_range_str.c_str()) ;

  ExpPlot.SetTitle(tmp) ;
    
  TCanvas c_exp("c_exp_e", "ExpVSe", 1500, 1000) ;
  c_exp.SetGrid();

  gPad->SetGrid() ;
  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.035);
  gStyle->SetTitleAlign(23) ;
  c_exp.Update() ;

  ExpPlot.SetNpx(5000);

  ExpPlot.Draw() ;

  c_exp.SaveAs((wrk_dir + f_name +".pdf").Str().c_str()) ;

  Z_LOG_INFO("Exposure vs energy plot saved into: '"+f_name.Str()+".pdf'.") ;

}

//--------------------------------------------------------------
void Satellite::PlotExposureTime(const Zaki::String::Directory& f_name, double energy)
{
  if(!set_exposure_eval_flag)
  {
    Z_LOG_ERROR("Plot failed because exposure hasn't been evaluated!");
    exit(EXIT_FAILURE) ;
    return ;
  }

  Z_LOG_INFO("Plotting the exposure vs. time...") ;

  TH1F h = *(TH1F*)(exp_hist.Clone(("exposure_e_"+GetName()).c_str()));

  double scale_factor = 24*3600*ExpTimeFrac(energy);
  scale_factor       *= Acceptance(energy) / GetConeFOV("sr");
  h.Scale(scale_factor) ;

  // h.GetXaxis()->SetTitle("Day");

  std::string date_range_str ;
  date_range_str = GetTimeDuration().start.StrForm("#Y-#M-#D")
                  + "   to   "
                  + GetTimeDuration().end.StrForm("#Y-#M-#D");
  char tmp[200] ;
  sprintf(tmp, "#splitline{Daily Exposure(E = %.1f  GeV) vs. Time}{%s}; t (day); Exposure (m^2.s)",
          energy, date_range_str.c_str()) ;

  h.SetTitle(tmp);
  
  TCanvas c_exp("c_exp", "Exp", 2500, 1000) ;
  c_exp.SetGrid();

  gPad->SetGrid() ;
  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.035);
  gStyle->SetTitleAlign(23) ;
  c_exp.Update() ;

  h.GetXaxis()->SetTimeDisplay(1) ;
  TDatime da(JAN1st1995.yr,JAN1st1995.mo,JAN1st1995.day,
             JAN1st1995.cl.h,JAN1st1995.cl.min,JAN1st1995.cl.sec);
  h.GetXaxis()->SetTimeOffset(da.Convert()) ;
  h.GetXaxis()->SetTimeFormat("#splitline{%Y}{%b/%d}") ;
  h.GetXaxis()->SetLabelSize(0.03) ;
  h.GetXaxis()->SetLabelOffset(0.015) ;

  h.GetXaxis()->SetTitleOffset(1.5) ; 

  h.Draw("HIST") ;

  char out_file_char[150] ;    
  sprintf(out_file_char, "_%.0f.pdf", energy) ;

  c_exp.SaveAs((wrk_dir + f_name + std::string(out_file_char)).Str().c_str()) ;

  Z_LOG_INFO("Exposure plot saved into: '" + f_name.Str() + std::string(out_file_char) + "'.") ;
}

//--------------------------------------------------------------
void Satellite::PlotSigShape(const Zaki::String::Directory& f_name)
{
  // sig_shape_hist.GetXaxis()->SetTitle("Day");
  // sig_shape_hist.GetYaxis()->SetTitle("Normalized Exposure");

  // sig_shape_hist.SetTitle("Signal Shape") ;

  std::string date_range_str ;
  date_range_str = GetTimeDuration().start.StrForm("#Y-#M-#D")
                  + "   to   "
                  + GetTimeDuration().end.StrForm("#Y-#M-#D");
  char tmp[150] ;
  sprintf(tmp, "#splitline{Signal Shape vs. Time}{%s}; t (day); Normalized Exposure",
          date_range_str.c_str()) ;

  sig_shape_hist.SetTitle(tmp) ;

  TCanvas c_sig("c_sig", "AMS-02", 1000, 600) ;
  c_sig.SetGrid();

  gPad->SetGrid() ;
  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.035);
  gStyle->SetTitleAlign(23) ;
  c_sig.Update() ;

  sig_shape_hist.Draw("HIST") ;

  c_sig.SaveAs((wrk_dir + f_name + ".pdf").Str().c_str()) ;
}

//--------------------------------------------------------------
void Satellite::ResetSigShape()
{
  if(!set_exposure_eval_flag)
  {
    Z_LOG_ERROR(" ResetSigShape() failed because exposure hasn't been evaluated!");
    exit(EXIT_FAILURE) ;
    return ;
  }

  TH1F *h1 = (TH1F*)(exp_hist.Clone(("sig_shape_"+GetName()).c_str()));
  sig_shape_hist = *h1 ;
}

//--------------------------------------------------------------
void Satellite::NormalizeExposure() 
{ 
  if(!set_exposure_eval_flag)
  {
    Z_LOG_ERROR(" NormalizeExposure() failed because exposure hasn't been evaluated!");
    exit(EXIT_FAILURE) ;
    return ;
  }

  TH1F *h1 = (TH1F*)(exp_hist.Clone(("sig_shape_"+GetName()).c_str()));
  double scale_factor = 1.0 / h1->Integral();

  h1->Scale(scale_factor) ;

  sig_shape_hist = *h1 ;
}

//--------------------------------------------------------------
void Satellite::ImportData(const std::string& in_name, const Zaki::String::Directory& file_name, 
                    const std::vector<DMSolarSignal::Particle>& in_states)
{
  results.SetStates(in_states);
  results.Input(file_name)  ;
  results.SetName(in_name) ;
  set_data_flag = true ;

  if(set_wrk_dir_flag)
    results.SetWrkDir(wrk_dir) ;
}

//--------------------------------------------------------------
Data Satellite::GetData() const
{
  return results;
}

//--------------------------------------------------------------
double Satellite::GetNumScaling() const
{
  return NUM_SCALING;
}

//--------------------------------------------------------------
double Satellite::GetExpNorm() const
{
    if(!set_exposure_eval_flag)
  {
    Z_LOG_ERROR(" GetExpNorm failed because exposure hasn't been evaluated!");
    exit(EXIT_FAILURE) ;
    return -1;
  }

  return exp_hist.Integral() ; 
}

//--------------------------------------------------------------
void Satellite::TimeBin(int bin_period)
{
  ResetSigShape();
  sig_shape_hist.Rebin(bin_period);

  double scale_factor = 1.0 / sig_shape_hist.Integral();
  sig_shape_hist.Scale(scale_factor) ;
}

//--------------------------------------------------------------
TH1F Satellite::GetSigShape(double energy) const
{
  TH1F h = *(TH1F*)(sig_shape_hist.Clone(("exposure_e_"+GetName()).c_str()));

  double scale_factor = NUM_SCALING*24*3600*ExpTimeFrac(energy);
  scale_factor       *= Acceptance(energy) / GetConeFOV("sr");
  scale_factor       *= GetExpNorm() ;

  h.Scale(scale_factor) ;

  return h;
}

//--------------------------------------------------------------
void Satellite::SetOrbitInclination(const double& o)
{
  set_orbit_inclination_flag = true;
  orbit_inclination = o;
}

//--------------------------------------------------------------
double Satellite::GetOrbitInclination() const
{
  if (set_orbit_inclination_flag)
    return orbit_inclination ;
  else
  {
    Z_LOG_ERROR("Orbit inclination hasn't been set!") ;
    return -1 ;
  }
}

//--------------------------------------------------------------
Zaki::Physics::DateInterval Satellite::GetTimeDuration() const 
{
  if (set_time_duration_flag)
    return t_duration ;
  else 
  {
    Z_LOG_ERROR("Time duration hasn't been set!") ;
    return {{2000, 1, 1, {0, 0, 0}}, {2000, 1, 1, {0, 0, 0}}} ;
  }
}

//--------------------------------------------------------------
Zaki::Physics::DateInterval Satellite::GetActualTimeDuration() const 
{
  // if (set_time_duration_flag)
    return t_duration_actual ;
  // else 
  // {
    // Z_LOG_ERROR("Time duration hasn't been set!") ;
    // return {{2000, 1, 1, {0, 0, 0}}, {2000, 1, 1, {0, 0, 0}}} ;
  // }
}

//--------------------------------------------------------------
void Satellite::SetTimeDuration(const Zaki::Physics::DateInterval& in_dur)
{
  set_time_duration_flag = true ;
  t_duration = in_dur;

  results.ScaleData(t_duration.duration() / GetActualTimeDuration().duration()) ;
}

//--------------------------------------------------------------
void Satellite::SetActualTimeDuration(const Zaki::Physics::DateInterval& in_dur)
{
  // set_time_duration_flag = true ;
  t_duration_actual = in_dur;
}

//--------------------------------------------------------------
// void Satellite::SetEclipticSolarLongRate(double Gamma)
// {
//   set_ecliptic_solar_long_rate_flag = true;
//   ecliptic_solar_long_rate = Gamma ;
// }

//--------------------------------------------------------------
// Solar ecliptic longitude
// double Satellite::GetEclipticSolarLong(const double& t) const
// {
//   // (0.9856 Degree)/day (t + t0);
//   // return SOLAR_PREC_RATE * (t + std::get<0>(GetTimeDuration())) ;
//   return SOLAR_PREC_RATE * t ;
// }

//--------------------------------------------------------------
void Satellite::SetEccentricity(const double& ec)  
{
  set_eccentricity_flag    = true ;
  eccentricity = ec ;
}

//--------------------------------------------------------------
void Satellite::SetOrbitHeight(const double& h)   // km
{
  set_orbit_height_flag    = true ;
  orbit_height = h  ;

  // avoiding circular definitions
  if(!set_orbit_radius_flag)  
    SetOrbitRadius( EARTH_RADIUS + h) ;

  // updating the radius if it's not:  r = h + R_earth
  else if (orbit_radius != EARTH_RADIUS + orbit_height)
  {
    SetOrbitRadius( EARTH_RADIUS + h) ;
  }
}

//--------------------------------------------------------------
void Satellite::SetOrbitRadius(const double& r)   // km
{
  set_orbit_radius_flag    = true ;
  orbit_radius = r ;

  // avoiding circular definitions
  if(!set_orbit_height_flag)  
    SetOrbitHeight(r - EARTH_RADIUS) ;
  // updating the height if it's not:  h = r - R_earth
  else if (orbit_height != orbit_radius - EARTH_RADIUS)
  {
    SetOrbitHeight(r - EARTH_RADIUS) ;
  }
}

//--------------------------------------------------------------
void Satellite::SetOrbitPeriod(const double& t)   // minutes
{
  set_orbit_period_flag    = true ;
  orbit_period = t ;

  // avoiding circular definitions 
  if(!set_orbit_ang_freq_flag)  
    SetOrbitAngFreq(2*M_PI/t) ; // (radian / minutes)
  // updating the orbit frequency if it's not:  2*M_PI/t
  else if (orbit_ang_freq != 2*M_PI/t)
  {
    SetOrbitAngFreq(2*M_PI/t) ;
  }
}

//--------------------------------------------------------------
void Satellite::SetOrbitAngFreq(const double& f) // (radian / minutes)
{
  set_orbit_ang_freq_flag  = true ;
  orbit_ang_freq = f  ;

  // avoiding circular definitions
  if(!set_orbit_period_flag)   
    SetOrbitPeriod(2*M_PI/f)  ; // minutes
  // updating the orbit period if it's not:  2*M_PI/f
  else if (orbit_period != 2*M_PI/f)
  {
    SetOrbitPeriod(2*M_PI/f) ;
  }
}

//--------------------------------------------------------------
double Satellite::GetEccentricity()  const
{
  if (set_eccentricity_flag)
    return eccentricity ;
  else 
    return -1 ;
}

//--------------------------------------------------------------
double Satellite::GetOrbitHeight()  const 
{
  if (set_orbit_height_flag)
    return orbit_height ;
  else 
    return -1 ;
}

//--------------------------------------------------------------
double Satellite::GetOrbitRadius()  const  
{
  if (set_orbit_radius_flag)
    return orbit_radius ;
  else 
    return -1 ;
}

//--------------------------------------------------------------
double Satellite::GetOrbitPeriod() const  
{
  if (set_orbit_period_flag)
    return orbit_period ;
  else 
    return -1 ;
}

//--------------------------------------------------------------
double Satellite::GetOrbitAngFreq() const 
{
  if (set_orbit_ang_freq_flag)
    return orbit_ang_freq ;
  else 
    return -1 ;
}

//--------------------------------------------------------------
// Forces  precession rate to match (360/365.25)
// to have a Sun synchronous orbit 
void Satellite::ForceSunSync()
{
  FindPrecessionRate() ;

  // Checking if the precession rate is more than 2 degrees
  // different from the Sun's precession rate.
  if (abs(precession_rate - 360.0/365.25) > 2)
  {
    Z_LOG_ERROR("The orbit is far from being synchronous, forcing failed.") ;
    Z_LOG_WARNING("The value from the evaluated precession rate will be used instead.") ;
    return ;
  }

  precession_rate = 360.0 / 365.25 ;
  forced_sun_sync_flag = true ;
  found_precession_rate_flag = false;
}

//--------------------------------------------------------------
// Precession of the orbital node:
// if it matches the apparent motion of the Sun (+1 degrees/day)
// it results in the Sun-synchronous orbit, note that inclination
// must be greater than 90, so that Cos[i] <0 to get a positive inclination.
void Satellite::FindPrecessionRate()
{ 
  //..........................................................
  //    Safety Checks
  int failure = 0 ;
  if (!set_orbit_radius_flag)
  {
    Z_LOG_ERROR("Orbit radius/height not set!") ;
    failure++ ;
  }
  if (!set_eccentricity_flag)
  {
    Z_LOG_ERROR("Orbit eccentricity not set!")  ;
    failure++ ;
  }
  if (!set_orbit_inclination_flag)
  {
    Z_LOG_ERROR("Orbit inclination not set!")   ;
    failure++ ;
  }
  if (!set_orbit_ang_freq_flag)
  {
    Z_LOG_ERROR("Orbit angular frequency or period not set!")   ;
    failure++ ;
  }

  if (failure != 0)
  {
    Z_LOG_ERROR("Missing a total of '"+ std::to_string(failure) + "' parameters.") ;
    return ;
  }
  //..........................................................

  double Re_sqr   = pow(EARTH_RADIUS, 2) ;
  double R_sqr    = pow(orbit_radius, 2) ;
  double ecc_sqr  = pow(eccentricity, 2) ;
  double cos_i    = cos(orbit_inclination*DEG_2_RAD)  ;


  precession_rate = -(3.0/2.0)*Re_sqr*EARTH_J2*orbit_ang_freq*cos_i / (R_sqr* pow(1 - ecc_sqr,2) ) ;

  // Converting from (rad/min) to (deg/day)
  precession_rate *= (RAD_2_DEG/MIN_2_DAY) ; 

  found_precession_rate_flag = true ;
}

//--------------------------------------------------------------
// in (deg / day)
double Satellite::GetPrecessionRate() const
{
  if (!found_precession_rate_flag && !forced_sun_sync_flag)
  {
    Z_LOG_ERROR("Precession rate hasn't been found yet, use 'FindPrecessionRate()' method first!") ;
    return -99999 ;
  }
  return precession_rate ;
}

//--------------------------------------------------------------
// RAAN at vernal equinox
void Satellite::SetInitRightAscNode(const double& right_node)
{
  init_right_asc_node = right_node ;
  set_right_asc_node_flag = true ;
}

//--------------------------------------------------------------
// RAAN at vernal equinox
double Satellite::GetInitRightAscNode() const
{
  if (set_right_asc_node_flag)
    return init_right_asc_node;
  else 
    return -1;
}

//--------------------------------------------------------------
// Local time descending node
// 't' in UnixTDays
Zaki::Physics::Clock Satellite::GetLocalTDescNode(double t) const  
{
  return GetLocalTAscNode(t) + 12 ;
}

//--------------------------------------------------------------
// Local time ascending node
// 't' in UnixTDays
Zaki::Physics::Clock Satellite::GetLocalTAscNode(double t) const   
{
  Zaki::Physics::Clock zc(12 + 24*(GetRAAN(t) -
                              GetSunPos(t).alpha()) / 360) ;

  // zc = {12 + 24*(GetSatPos(t).ToEcliptic(t).lambda -
  //   GetSunEclipticLong(t)) / 360};

  return zc ;
}

//--------------------------------------------------------------
// RAAN 
double Satellite::GetRAAN(double in_Jday) const 
{
  double tmp_time_dur = in_Jday - GetTimeDuration().start.UnixTDay() ;
  double tmp_RAAN = GetInitRightAscNode() + GetPrecessionRate()*tmp_time_dur ;

  return tmp_RAAN - 360*floor( tmp_RAAN / 360 ) ;
}

//--------------------------------------------------------------
// Initial theta (True Anomaly) in degrees
void Satellite::SetInitTrueAnomaly(const double& theta_0)
{
  init_true_anomaly = theta_0 ;
  set_init_true_anomaly_flag = true ;
}

//--------------------------------------------------------------
// True Anomaly Initial True Anomaly
double Satellite::GetInitTrueAnomaly() const
{
  if (set_init_true_anomaly_flag)
    return init_true_anomaly;
  else 
    return -1;
}

//--------------------------------------------------------------
// True Anomaly (deg)
double Satellite::GetTrueAnomaly(double in_date) const // 't' in UnixTDays
{
  double tmp_time_dur = in_date - GetTimeDuration().start.UnixTDay() ;
  double tmp_theta = GetInitTrueAnomaly() +
                     GetOrbitAngFreq()*tmp_time_dur*(RAD_2_DEG/MIN_2_DAY);

  return tmp_theta - 360*floor( tmp_theta / 360 ) ;
}

//--------------------------------------------------------------
// Solar Beta Angle ( t in UnixTDays)
double Satellite::SolarBeta(double t) const
{
  double Gamma = GetSunEclipticLong(t)*DEG_2_RAD    ;
  double i     = GetOrbitInclination()*DEG_2_RAD    ;
  double eps   = Zaki::Physics::EclipticObliquity(t)*DEG_2_RAD ;
  double omega = GetRAAN(t)*DEG_2_RAD               ;

  double beta = ( cos(Gamma)*sin(omega)*sin(i)
                  - sin(Gamma)*cos(eps)*cos(omega)*sin(i)
                  + sin(Gamma)*sin(eps)*cos(i) ) ;

  beta = asin(beta)*RAD_2_DEG ;

  return beta;
}

//--------------------------------------------------------------
// Wrapper for making plots using root
double Satellite::BetaFuncWrapper(double *x, double *par)
{
  double tmp = (*x/86400.0) + JAN1st1995.UnixTDay() ;
  return SolarBeta(tmp);
}

//--------------------------------------------------------------
// 
// double Satellite::OrbitFunc(Zaki::Physics::DateInterval in_interval,
//    size_t in_num)
// {
//   std::vector<double> x, y, z ;
//   x.reserve(in_num+1) ; y.reserve(in_num+1) ; z.reserve(in_num+1) ;
//   double delta_t = in_interval.duration() / in_num ;

//   for (size_t i = 0; i < in_num+1; i++)
//   {
//     x.push_back(GetSatPos(in_interval.start.UnixTDay() + i*delta_t).xyz.x()) ;
//     y.push_back(GetSatPos(in_interval.start.UnixTDay() + i*delta_t).xyz.y()) ;
//     z.push_back(GetSatPos(in_interval.start.UnixTDay() + i*delta_t).xyz.z()) ;
//   }
  
//   return GetSatPos(1).alpha();
// }

//--------------------------------------------------------------
void Satellite::PlotOrbit(const Zaki::Physics::Date& t_1,
  const Zaki::Physics::Date& t_2, const Zaki::String::Directory& f_name)
{
  size_t in_num = 150 ;

  std::vector<double> x, y, z ;
  x.reserve(in_num+1) ; y.reserve(in_num+1) ; z.reserve(in_num+1) ;
  double delta_t = (t_2 - t_1) / in_num ;

  for (size_t i = 0; i < in_num+1; i++)
  {
    x.push_back(GetSatPos(t_1.UnixTDay() + i*delta_t).XYZ().x()) ;
    y.push_back(GetSatPos(t_1.UnixTDay() + i*delta_t).XYZ().y()) ;
    z.push_back(GetSatPos(t_1.UnixTDay() + i*delta_t).XYZ().z()) ;
  }

  TCanvas c_orbit("c_orbit","Orbit",0,0,800,800);

  // ...............................................................
  // Origin
  TGraph2D origin(2) ; origin.SetPoint(0, 0, 0, 0) ;  origin.SetPoint(0, 0, 0, 0.5) ;
  origin.SetMarkerSize(2) ; origin.SetMarkerStyle(20) ; origin.SetName("origin") ;
  // ...............................................................

  TGraph2D orbit_graph(static_cast<Int_t>(in_num), &x[0], &y[0], &z[0]) ;
  orbit_graph.SetName("Orbit") ;
  orbit_graph.SetMarkerStyle(20);

  orbit_graph.GetXaxis()->SetTitle("I");
  orbit_graph.GetXaxis()->SetLabelSize(0.02) ;
  orbit_graph.GetXaxis()->CenterTitle() ;
  orbit_graph.GetXaxis()->SetTitleOffset(2.0) ;
  orbit_graph.GetYaxis()->SetTitle("J");
  orbit_graph.GetYaxis()->SetLabelSize(0.02) ;
  orbit_graph.GetYaxis()->CenterTitle() ;
  orbit_graph.GetYaxis()->SetTitleOffset(2.5) ;
  orbit_graph.GetZaxis()->SetTitle("K");
  orbit_graph.GetZaxis()->SetLabelSize(0.02) ;
  orbit_graph.GetZaxis()->CenterTitle() ;
  orbit_graph.GetZaxis()->SetTitleOffset(1.5) ;

  char tmp_title[200] ;
  sprintf(tmp_title, "%s Orbit [%.1f, %.1f]", GetName().c_str(), t_1.UnixTDay(), t_2.UnixTDay()) ;
  orbit_graph.SetTitle(tmp_title) ;

  orbit_graph.SetLineColorAlpha(kBlue, 0.5) ;

  // gStyle->SetPalette(1);
  // orbit_graph.Draw("pcol");
  // earth_sphere.Draw("p0") ;
  orbit_graph.Draw("line");
  origin.Draw("P same") ;


  c_orbit.SaveAs((wrk_dir + f_name + ".pdf").Str().c_str()) ;

  Z_LOG_INFO("Satellite orbit plot saved into: '"+f_name.Str()+"'.") ;
}
//--------------------------------------------------------------
// Plotting Solar Beta Angle from t_1 to t_2
void Satellite::PlotSolarBeta(const Zaki::Physics::Date& t_1,
  const Zaki::Physics::Date& t_2, const Zaki::String::Directory& f_name)
{
  TCanvas c_beta("c_beta", "Beta Angle", 2500, 1400) ;
  c_beta.SetGrid();
  gPad->SetGrid() ;
  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.035);
  gStyle->SetTitleAlign(23) ;
  c_beta.Update() ;
  // TF1("f",fptr,&MyFunction::Evaluate,0,1,npar,"MyFunction","Evaluate"); 
  TF1 BetaPlot("Beta", this, &Satellite::BetaFuncWrapper, t_1.RootT95(), t_2.RootT95(), 0);

  std::string date_range_str ;
  date_range_str = t_1.StrForm("#Y-#M-#D")
                  + "   to   "
                  + t_2.StrForm("#Y-#M-#D");

  char tmp[150] ;
  sprintf(tmp, "#splitline{Solar Beta Angle vs. Time}{%s}; t (day); #beta(t)",
          date_range_str.c_str()) ;

  BetaPlot.SetTitle(tmp) ;
  BetaPlot.SetNpx(5500);

  
  Zaki::Math::Range<double> tmp_y = { 
      BetaPlot.GetMinimum(t_1.RootT95(), t_2.RootT95()) - 5,
      BetaPlot.GetMaximum(t_1.RootT95(), t_2.RootT95()) + 5 };

  BetaPlot.SetMaximum(tmp_y.max);
  BetaPlot.SetMinimum(tmp_y.min); 

  BetaPlot.GetXaxis()->SetTimeDisplay(1) ;
  TDatime da(JAN1st1995.yr,JAN1st1995.mo,JAN1st1995.day,
             JAN1st1995.cl.h,JAN1st1995.cl.min,JAN1st1995.cl.sec);
  BetaPlot.GetXaxis()->SetTimeOffset(da.Convert()) ;
  BetaPlot.GetXaxis()->SetTimeFormat("#splitline{%Y}{%b/%d}") ;
  BetaPlot.GetXaxis()->SetLabelSize(0.028) ;
  BetaPlot.GetXaxis()->SetLabelOffset(0.02) ;

  BetaPlot.GetXaxis()->SetTitleOffset(1.48) ;

  BetaPlot.Draw() ;

  c_beta.SaveAs((wrk_dir + f_name + ".pdf").Str().c_str()) ;

  Z_LOG_INFO("Solar beta angle plot saved into: '"+f_name.Str()+"'.") ;
}

//--------------------------------------------------------------
void Satellite::SetYPR(const double& y, const double& p, const double& r) 
{
  ypr.Yaw = y; ypr.Pitch = p ; ypr.Roll = r;
  set_YPR_flag = true;

  // The normal to the detector in the satellite local frame
  // which is non-inertial:
  // Yaw   ---->  Z'
  // Pitch ---->  Y'
  // Roll  ---->  X'
  normal.SetCoordinates(0, 0, -1) ;
  ROOT::Math::RotationZYX r_local(y*DEG_2_RAD, p*DEG_2_RAD, r*DEG_2_RAD);
  normal = r_local*normal ;
}

//--------------------------------------------------------------
void Satellite::SetYPR(const YPR& in_ypr) 
{
  ypr = in_ypr ;
  set_YPR_flag = true;

  // The normal to the detector in the satellite local frame
  // which is non-inertial:
  // Yaw   ---->  Z'
  // Pitch ---->  Y'
  // Roll  ---->  X'
  // (X', Y', Z') = (0, 0, 0) is the position of the satellite
  normal.SetCoordinates(0, 0, -1) ;
  ROOT::Math::RotationZYX r_local(ypr.Rad().Yaw, 
                                  ypr.Rad().Pitch,
                                  ypr.Rad().Roll);
  normal = r_local*normal ;
}

//--------------------------------------------------------------
// The normal to the detector in the satellite local frame
// which is non-inertial:
// Yaw   ---->  Z'
// Pitch ---->  Y'
// Roll  ---->  X'
// (X', Y', Z') = (0, 0, 0) is the position of the satellite
ROOT::Math::XYZVector Satellite::GetNormLocal(double in_JDate) const
{
  return normal;
}

//--------------------------------------------------------------
// The normal to the detector in the satellite GEI frame
ROOT::Math::XYZVector Satellite::GetNormGEI(double in_JDate) const
{
  // Z' is in the nadir direction, i.e. '-R_Sat' in the GEI coordinate
  ROOT::Math::XYZVector Loc_X = GetSatTan(in_JDate).XYZ()  ;
  ROOT::Math::XYZVector Loc_Z = -GetSatPos(in_JDate).XYZ() /
                                 sqrt(GetSatPos(in_JDate).XYZ().mag2()) ;
  ROOT::Math::XYZVector Loc_Y = Loc_Z.Cross(Loc_X) ;
  Loc_Y = Loc_Y / sqrt(Loc_Y.mag2());

  ROOT::Math::XYZVector GEI_X = {1, 0, 0} ;
  ROOT::Math::XYZVector GEI_Y = {0, 1, 0} ;
  ROOT::Math::XYZVector GEI_Z = {0, 0, 1} ;

  ROOT::Math::XYZVector GEI_normal = normal ;

  // Rotation matrix R for
  // GEI --> R --> Local(primed)
  // i.e. v'_i = R_i^j*v_j
  // so R_i^j = < e_i' | e^j >
  std::array<double, 9> rotAr ;  
  rotAr[0] = Loc_X.Dot(GEI_X); rotAr[1] = Loc_X.Dot(GEI_Y); rotAr[2] = Loc_X.Dot(GEI_Z); 
  rotAr[3] = Loc_Y.Dot(GEI_X); rotAr[4] = Loc_Y.Dot(GEI_Y); rotAr[5] = Loc_Y.Dot(GEI_Z); 
  rotAr[6] = Loc_Z.Dot(GEI_X); rotAr[7] = Loc_Z.Dot(GEI_Y); rotAr[8] = Loc_Z.Dot(GEI_Z);

  ROOT::Math::Rotation3D rot(rotAr.begin(), rotAr.end()) ;
  // .................................
  //      Euler angles
  //
  //    GEI  ---> EulerRotation ---> Local(shifted)
  //    Local(shifted)  ---> Invert(EulerRotation) ---> GEI
  //
  // double EU_theta, EU_phi, EU_psi;
  
  // if ( r22 < +1)
  // {
  //   if ( r22 > −1)
  //   {
  //     EU_theta = acos ( r22 ) ;
  //     EU_psi = atan2 ( r02 , r12 ) ;
  //     EU_phi = atan2 ( r20 , -r21 ) ;
  //   }
  //   else // r22 = −1
  //   {
  //     // Not a unique solution: psi - phi  = atan2 (−r01 , r00 )
  //     EU_theta = M_PI ;
  //     EU_psi = atan2 (−r01 , r00 ) ;
  //     EU_phi = 0;
  //   }
  // }
  // else // r22 = +1
  // {
  //   // Not a unique solution: phi + psi = atan2 (r01 , r00 )
  //   EU_theta = 0;
  //   EU_psi = atan2 (r01 , r00 ) ;
  //   EU_phi = 0;
  // }

  // Shifting from local frame to the GEI origin:
  // GEI_normal = GEI_normal + tmp_sat_pos ;
  // Finding the line of nodes: N = (Z' X Z), where Z' is in the 
  // nadir direction, i.e. '-R_Sat' in the GEI coordinate
  // ROOT::Math::XYZVector l_nodes = GEI_Z.Cross(Loc_Z) ;

  // std::cout << "GEI_Z.Cross(Loc_Z): " << l_nodes << "\n" ;
  // std::cout << "Loc_Z.Cross(GEI_Z): " << Loc_Z.Cross(GEI_Z) << "\n" ;
  //
  // Theta : between Z & Z'
  // double EU_theta = -acos(Loc_Z.Dot(GEI_Z) /
  //                        sqrt(GEI_Z.mag2()) / sqrt(Loc_Z.mag2()) );

  // std::cout << "EU_theta: " << EU_theta << "\n" ;
  // //
  // // phi : between X & N
  // double EU_phi = acos(l_nodes.Dot(GEI_X) /
  //                       sqrt(GEI_X.mag2()) / sqrt(l_nodes.mag2()) );
  // std::cout << "EU_phi: " << EU_phi << "\n" ;
  // //
  // // psi : between X' & N
  // double EU_psi = acos(l_nodes.Dot(GetSatTan(in_JDate).xyz) /
  //                       sqrt(GetSatTan(in_JDate).xyz.mag2()) /
  //                       sqrt(l_nodes.mag2()) );
  // std::cout << "EU_psi: " << EU_psi << "\n" ;
  //
  // EulerAngles ( phi,  theta,  psi) :
  // ROOT::Math::EulerAngles EU_ang( EU_phi, EU_theta, EU_psi );
  // .................................

  // .................................
  // Rotating the local frame back onto the GEI:
  //
  //    GEI  ---> EulerRotation ---> Local(shifted)
  //    Local(shifted)  ---> Invert(EulerRotation) ---> GEI
  //
  // EU_ang.Invert() ;
  // ROOT::Math::Rotation3D EU_Rot_3D(EU_ang) ;
  // GEI_normal = EU_Rot_3D*GEI_normal ;
  // .................................

  // Shifting back from GEI origin to local frame:
  // GEI_normal = GEI_normal - tmp_sat_pos ;
  
  GEI_normal = rot.Inverse()*GEI_normal ;
  return GEI_normal / sqrt(GEI_normal.mag2()) ;
}

//--------------------------------------------------------------
YPR Satellite::GetYPR() const
{
  if(set_YPR_flag)
    return ypr ;
  else
    return {-1, -1, -1} ;
}
//--------------------------------------------------------------
// The vector from the Sun to the center of the Satellite orbit
//  in the inertial orbit coordinate. 
// ROOT::Math::XYZVector Satellite::GetRSun2C(double t) const
// {
//   ROOT::Math::XYZVector v(0,  EARTH_2_SUN*cos(Beta(t)*DEG_2_RAD),
//                               -EARTH_2_SUN*sin(Beta(t)*DEG_2_RAD)) ;

//   return v;
// }

//--------------------------------------------------------------
// The vector from the Earth's center to the satellite 
// in the inertial orbit coordinate. 
// ROOT::Math::XYZVector Satellite::GetRC2Sat(double t) const
// {
//   ROOT::Math::XYZVector v ;
//   double ra = GetOrbitHeight() ;
//   double thet = GetTrueAnomaly(t)*DEG_2_RAD ;

//   v.SetCoordinates(ra*sin(thet), - ra*cos(thet), 0) ;

//   return v;
// }

//--------------------------------------------------------------
// The vector from the Sun to the Satellite in the inertial
//  orbit coordinate. 
Zaki::Physics::GEICoord  Satellite::GetRSun2Sat(double t) const
{
  return -GetSunPos(t) + GetSatPos(t);
}

//--------------------------------------------------------------
// The projection of the vector from the Sun on to the Satellite
double Satellite::GetSun2SatProj(double t) const
{
  return GetRSun2Sat(t).XYZ().Dot(GetNormGEI(t)) / sqrt(GetRSun2Sat(t).XYZ().mag2()) ;
}

//--------------------------------------------------------------
// Plotting Solar projection "cos(r_Sun.n_ams)" from t_1 to t_2
void Satellite::PlotSunProj(const Zaki::Physics::Date& t_1,
const Zaki::Physics::Date& t_2, const Zaki::String::Directory& f_name)
{
  // ............................................
  // Exporting the data points
  std::vector<double> sunshine_vec ;
  sunshine_vec.reserve(4*(t_2.JD() - t_1.JD())) ;
  for(double i=t_1.RootT95() ; i <= t_2.RootT95() ; i += 3600*6 )
  {
    sunshine_vec.push_back(SunshineFuncWrapper(&i, nullptr)) ;
  }
  Zaki::File::VecSaver my_saver(wrk_dir + f_name + ".dat", Zaki::File::FileMode::Write) ;
  my_saver.Export1D(sunshine_vec) ;
  // ............................................

  TCanvas c_sunshine("c_sunshine", "Sunshine Projection", 2500, 1200) ;
  c_sunshine.SetGrid();
  gPad->SetGrid() ;
  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.035);
  gStyle->SetTitleAlign(23) ;
  c_sunshine.Update() ;

  // TF1("f",fptr,&MyFunction::Evaluate,0,1,npar,"MyFunction","Evaluate"); 
  TF1 SunshinePlot("Sunshine", this, &Satellite::SunshineFuncWrapper, t_1.RootT95(), t_2.RootT95(), 0);

  std::string date_range_str ;
  date_range_str = GetTimeDuration().start.StrForm("#Y-#M-#D")
                  + "   to   "
                  + GetTimeDuration().end.StrForm("#Y-#M-#D");
  char tmp[150] ;
  sprintf(tmp, "#splitline{Sunshine on the Satellite vs. Time}{%s}; t (day); #hat{r}_{Sun} #bullet #hat{n}",
          date_range_str.c_str()) ;

  SunshinePlot.SetTitle(tmp) ;

  SunshinePlot.SetNpx(35000);
  SunshinePlot.SetMaximum(1);
  SunshinePlot.SetMinimum(-1);
  SunshinePlot.SetLineColor(kOrange) ;
  SunshinePlot.SetLineWidth(1) ;
  // SunshinePlot.SetLineStyle(7);
  // SunshinePlot.SetLineColorAlpha(0, 1) ;

  SunshinePlot.GetXaxis()->SetTimeDisplay(1) ;
  TDatime da(JAN1st1995.yr,JAN1st1995.mo,JAN1st1995.day,
             JAN1st1995.cl.h,JAN1st1995.cl.min,JAN1st1995.cl.sec);
  SunshinePlot.GetXaxis()->SetTimeOffset(da.Convert()) ;
  SunshinePlot.GetXaxis()->SetTimeFormat("#splitline{%Y}{%b/%d}") ;
  SunshinePlot.GetXaxis()->SetLabelSize(0.03) ;
  SunshinePlot.GetXaxis()->SetLabelOffset(0.015) ;

  SunshinePlot.GetXaxis()->SetTitleOffset(1.5) ;

  SunshinePlot.Draw() ;

  TLine line(t_1.RootT95(), -cos(GetConeFOV()*DEG_2_RAD),
             t_2.RootT95(), -cos(GetConeFOV()*DEG_2_RAD)) ;
  line.SetLineColor(kPink) ;
  line.SetLineStyle(2) ;
  line.Draw("SAME") ;

  c_sunshine.SaveAs((wrk_dir + f_name + ".pdf").Str().c_str()) ;

  Z_LOG_INFO("Solar projection on the satellite saved into: '"+f_name.Str()+"'.") ;
}

//--------------------------------------------------------------
// Wrapper for making plots using root
double Satellite::SunshineFuncWrapper(double *x, double *par)
{
  // converting from Root95 time to UnixTDay
  double tmp = (*x/86400.0) + JAN1st1995.UnixTDay() ;

  if ( IsInShadow(tmp) )
    return 0 ;
  else
    return GetSun2SatProj(tmp);
}

//--------------------------------------------------------------
// Input should be UnixTDays
double Satellite::ExposureIntegrand(double t)
{
  // (Sun -> Satellite) vector projection normal to the satellite's surface
  double tmp_proj = GetSun2SatProj(t) ;

  // Cosine of the field-of-view
  double cos_fov  = cos(GetConeFOV()*DEG_2_RAD);

  // Checking if within FOV & not in Earth's shadow
  if (tmp_proj < -cos_fov &&  !IsInShadow(t) )
    return -tmp_proj ;
  else
    return 0 ;
}

//--------------------------------------------------------------
bool Satellite::SunInFOV(const double& t)
{
  return GetSun2SatProj(t) < -cos(GetConeFOV()*DEG_2_RAD);
}

//--------------------------------------------------------------
// Evaluating the exposure between t_1 and t_2
void Satellite::EvaluateExposure(const Zaki::Physics::Date& t_1, const Zaki::Physics::Date& t_2) 
{
  TStopwatch timer;

  char tmp_char[200] ;
  sprintf(tmp_char, "Evaluating the exposure for %.1f days ...", t_2.UnixTDay()-t_1.UnixTDay()) ;
  Z_LOG_INFO(tmp_char) ;

  gsl_integration_workspace *w = gsl_integration_workspace_alloc(500*(size_t)(t_2.UnixTDay()-t_1.UnixTDay()+1));
  double err;

  Zaki::Math::GSLFuncWrapper<Satellite, double (Satellite::*)(double)> 
    Fp(this, &Satellite::ExposureIntegrand);     

  // ROOT::Math::WrappedMemFunction<Satellite, double ( Satellite::* ) (double)>
  //  f(*this, &Satellite::ExposureIntegrand);

  gsl_function F = static_cast<gsl_function> (Fp) ; 

  timer.Start();

  // ROOT::Math::Integrator ig(f);
  double s1 = 0, one_day = 0;
  std::vector<double> tmp_exp_set ;
  tmp_exp_set.reserve(static_cast<size_t>(t_2.UnixTDay()-t_1.UnixTDay())) ;

  // Finding the daily exposure and saving it in 'tmp_exp_set'
  for(size_t i = 0 ; i < t_2.UnixTDay()-t_1.UnixTDay() ; ++i)
  {
    gsl_integration_qag(&F, t_1.UnixTDay() + i, t_1.UnixTDay() + i + 1, 1e-3, 1e-3, 500, 1, w, &one_day, &err);
    // one_day = ig.Integral(t_1.UnixTDay() + i, t_1.UnixTDay() + i + 1) ;

    // weighing the events by ( 1 AU / R_sun_2_sat)^2 
    // Note that there should be a (1/AU)^2 in P_dec from the model
    // such that after cancellation we have included a time-dependent
    //  (1 / R_sun_2_sat )^2 factor.
    one_day *= AU_2_KM*AU_2_KM / pow(GetRSun2Sat(t_1.UnixTDay()+i).r(), 2) ; 

    s1 += one_day ;
    tmp_exp_set.push_back(one_day) ;
  }
  
  gsl_integration_workspace_free(w);

  Zaki::Math::Range<double> t_range  = { 0 , t_2.UnixTDay()-t_1.UnixTDay()} ;
  int bin_num   = static_cast<int> (t_2.UnixTDay()-t_1.UnixTDay()) ;

  TH1F tmp_exp_hist(("exp_"+GetName()).c_str(), 
                    ("Exposure_"+GetName()).c_str(), 
                    bin_num, t_range.min, t_range.max);

  // Filling the Exposure histogram
  for(size_t i=0; i<tmp_exp_set.size(); i++)
  {
    tmp_exp_hist.SetBinContent(i+1, tmp_exp_set[i])  ;
  }

  exp_hist = tmp_exp_hist;

  timer.Stop();
  
  double scale_exp = 24*3600*ExpTimeFrac(1.2e3);
  scale_exp       *= Acceptance(1.2e3) / GetConeFOV("sr");
  
  sprintf(tmp_char, "--> Integration took %f seconds.", timer.RealTime() ) ;
  Z_LOG_INFO(tmp_char) ;
  sprintf(tmp_char, "--> Integral result =  %f. \t Total Exposure (1.2 TeV): %.5e.", s1, s1*scale_exp ) ;
  Z_LOG_INFO(tmp_char) ;

  set_exposure_eval_flag = true ;
}

//--------------------------------------------------------------
void Satellite::ExportExposure(const Zaki::String::Directory& f_name, 
                               const Zaki::File::FileMode& mode) 
{
  if(!set_exposure_eval_flag)
  {
    Z_LOG_ERROR("Export failed because exposure hasn't been evaluated!");
    exit(EXIT_FAILURE) ;
    return ;
  }

  std::vector<double> exp_vec ;
  exp_vec.reserve(exp_hist.GetEntries()) ;

  for(size_t i=0 ; i < exp_hist.GetEntries() ; ++i)
    exp_vec.push_back(exp_hist.GetBinContent(i)) ;

  Zaki::File::VecSaver my_saver(wrk_dir + f_name, mode) ;
  my_saver.Export1D(exp_vec) ;

  Z_LOG_INFO("Exposure exported to '"+(wrk_dir +f_name).Str()+"'.") ;
}

//--------------------------------------------------------------
double Satellite::ExpTimeFrac(double energy) const
{
  return 1 ;
}

//--------------------------------------------------------------
double Satellite::Acceptance(double energy)  const
{
  return 1 ;
}

//--------------------------------------------------------------
double Satellite::GetMultiplicity(const std::vector<DMSolarSignal::Particle>& in_states) const
{
  if(!set_data_flag)
  {
    Z_LOG_ERROR("Data hasn't been imported yet!") ;
    return 0 ;
  }

  int counter = 0 ;
  for (auto s : in_states) 
  {
    if( Zaki::Vector::Exists(s, results.GetStates()) )
    {  
      counter++ ;
      // break;
    }
  }

  return counter;
}

//--------------------------------------------------------------
//==============================================================
