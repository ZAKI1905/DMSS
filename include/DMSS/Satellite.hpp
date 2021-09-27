#ifndef DMSS_Satellite_H
#define DMSS_Satellite_H

// Root
#include <Math/Vector3D.h>
#include <Math/GenVector/Rotation3D.h>
#include <Math/GenVector/EulerAngles.h>
#include <Math/GenVector/RotationX.h>
#include <Math/GenVector/RotationY.h>
#include <Math/GenVector/RotationZ.h>
#include <Math/GenVector/RotationZYX.h>

#include <Zaki/Physics/DateTime.hpp>
#include <Zaki/Physics/Constants.hpp>
#include <Zaki/Physics/Coordinate.hpp>

#include <Zaki/File/VecSaver.hpp>

// Local headers
#include "DMSS/Data.hpp"

//==============================================================
class Satellite : public Prog
{
  //--------------------------------------------------------------
  public:

    /// Constructor
    Satellite() ;

    /// Virtual destructor
    virtual ~Satellite() ;

    std::unique_ptr<Satellite> Clone() const
    {
      /// The derived class clones itself
      Satellite* tmp = IClone() ;
      
      /// We deep-copy the base-members containing
      /// pointer-members in this base class
      /// So no need to rewrite this in every copy
      /// constructor of the derived class
      tmp->exp_hist = *(static_cast<TH1F*> (tmp->exp_hist.Clone()));
      tmp->sig_shape_hist = *(static_cast<TH1F*> (tmp->sig_shape_hist.Clone()));
      
      /// returning a unique pointer to the cloned satellite
      return std::unique_ptr<Satellite>(tmp);
    }

    //.....................................
    // Setters
    //.....................................
    /// Sets the work directory
    void SetWrkDir(const Zaki::String::Directory&) override ;

    /// Sets the orbit inclination
    void SetOrbitInclination(const double&)           ;

    /// Sets the time period of data taking
    void SetTimeDuration(const Zaki::Physics::DateInterval&) ;
    void SetActualTimeDuration(const Zaki::Physics::DateInterval&) ;

    /// Filed of view (FOV) defined as the half-angle from the normal
    /// to the surface of the detector and in a conical shape.
    void SetConeFOV(const double&)             ;

    /// Sets the eccentricity of the orbit
    void SetEccentricity(const double&)        ;
    
    /// Sets the height of the orbit
    void SetOrbitHeight(const double&)         ;

    /// Sets the orbit of the orbit
    void SetOrbitRadius(const double&)         ;

    /// Sets the period of the orbit
    void SetOrbitPeriod(const double&)         ;

    /// Sets the angular frequency of the orbit
    void SetOrbitAngFreq(const double&)        ;

    /// Sets the right ascension of the ascending node at t_0
    void SetInitRightAscNode(const double&)    ;

    /// Sets true anomaly at t_0
    void SetInitTrueAnomaly(const double&)     ; 

    /// Sets orientation of the satellite on the orbit (Yaw, Pitch, Roll)
    void SetYPR(const double&, const double&, const double&) ; 

    /// Sets orientation of the satellite on the orbit (Yaw, Pitch, Roll)
    void SetYPR(const Zaki::Physics::YPR&) ;

    /// Imports the exposure data
    void ImportExposure(const Zaki::String::Directory&) ;

    /// Imports the satellite data
    void ImportData(const std::string&,
                    const Zaki::String::Directory&,
                    const std::vector<DMSolarSignal::Particle>&) ;

    //.....................................

    /// Checks whether the satelite is in the Earth's shadow
    bool IsInShadow(const double&); 

    /// Checks whether the Sun is in the Satellites's FOV
    bool SunInFOV(const double& ) ; 
    
    void TimeBin(int bin_period) ;

    /// Evaluates the exposure to the Sun given initial & final dates
    void EvaluateExposure(const Zaki::Physics::Date&, const Zaki::Physics::Date&) ;

    /// Exports the exposure into a file
    void ExportExposure(const Zaki::String::Directory&, 
          const Zaki::File::FileMode& = Zaki::File::FileMode::Write) ;
    
    /// The exposure integrand
    double ExposureIntegrand(double);

    /// Finds the precession rate
    void FindPrecessionRate();

    /// Forces  precession rate to match (360/365.25)
    /// to have a Sun synchronous orbit 
    void ForceSunSync();

    /// Normalizes the exposure histogram
    void NormalizeExposure();

    /// Resets the signal shape histogram
    void ResetSigShape();

    //.....................................
    // Getters
    //.....................................
    
    /// Sun position in GEI coordinates
    Zaki::Physics::GEICoord GetSunPos(double) const ;

    ///  The ecliptic longitude of the Sun
    double GetSunEclipticLong(double) const ;

    /// Satellite position in GEI coordinates, i.e. '-Yaw' direction
    Zaki::Physics::GEICoord GetSatPos(double) const ;

    /// Tangent to the orbit at satellite's position in GEI coordinates
    /// In the direction of motion, i.e. '+Roll' direction
    /// It's norm = 1
    Zaki::Physics::GEICoord GetSatTan(double) const ;

    /// Solar beta angle at time 't'
    double SolarBeta(double t)   const;
    
    /// Returns the daily exposure set
    std::vector<double> GetExposure() const ;

    /// Returns the signal shape histogram
    TH1F GetSigShape(double energy)   const ;

    /// returns the satellite's data
    Data GetData()                    const ;

    /// see SetKappa()
    double GetMultiplicity(const std::vector<DMSolarSignal::Particle>&) const ;
    
    // Satellite Parameters
    double GetConeFOV(std::string = "deg")  const ;
    double GetOrbitInclination()            const ;
    std::string GetOrbitType()              const ;
    Zaki::Physics::DateInterval GetTimeDuration()       const ;
    Zaki::Physics::DateInterval GetActualTimeDuration() const ;

    double GetEccentricity()              const ;
    double GetOrbitHeight()               const ;
    double GetOrbitRadius()               const ;
    double GetOrbitPeriod()               const ;
    double GetOrbitAngFreq()              const ;
    double GetPrecessionRate()            const ;
    double GetInitRightAscNode()          const ;

    /// Returns the local time descending node 
    Zaki::Physics::Clock GetLocalTDescNode(double)      const ;

    // Returns the local time ascending node
    Zaki::Physics::Clock GetLocalTAscNode(double)       const ; 
    double GetInitTrueAnomaly()                 const ;
    double GetNumScaling()                const ;

    /// Returns the normal to the detector's surface in GEI Coordinate
    ROOT::Math::XYZVector GetNormGEI(double)   const ;

    ///  Returns the normal to the detector's surface in Local coordinate (Y-P-R)
    ROOT::Math::XYZVector GetNormLocal(double) const ;

    Zaki::Physics::GEICoord GetRSun2Sat(double) const ;
    double GetSun2SatProj(double) const ;

    /// Returns the normalized exposure
    double GetExpNorm() const ;

    /// Returns the Right Ascention of the Ascending Node (RAAN) at time 't'   
    double GetRAAN(double)    const; 
    
    ///  Returns the true anomaly coordinate at time 't'  
    double GetTrueAnomaly(double t)  const;

    /// Returns the (Yaw, Pitch, Roll) of the satellite
    Zaki::Physics::YPR GetYPR() const ;
    //.....................................

    //.....................................
    // Plotters
    //.....................................
    /// Printing the satellite info
    void Print() const override;

    /// Plot the exposure over time & energy
    void PlotExposure(const Zaki::String::Directory&, double) ;

    /// Plot the exposure over time
    void PlotExposureTime(const Zaki::String::Directory&, double) ;

    /// Plot the total exposure over energy
    void PlotExposureEnergy(const Zaki::String::Directory&) ;

    /// Wrapper for: Satellite::PlotExposureEnergy
    double ExpEnergyFuncWrapper(double*, double*) ;

    /// Plot the normalized exposure
    void PlotSigShape(const Zaki::String::Directory&) ;

    /// Plots the solar beta angle from t_1 to t_2
    void PlotSolarBeta(const Zaki::Physics::Date&,
                       const Zaki::Physics::Date&, const Zaki::String::Directory&) ;

    /// Wrapper for Satellite::PlotSolarBeta
    double BetaFuncWrapper(double *x, double *par) ;

    /// Plots the projection of the Sun on the satellite
    void PlotSunProj(const Zaki::Physics::Date&,
                     const Zaki::Physics::Date&, const Zaki::String::Directory&) ;

    /// Wrapper for Satellite::PlotSunProj
    double SunshineFuncWrapper(double *, double *) ;

    /// Plots the orbit of the satellite
    void PlotOrbit(const Zaki::Physics::Date&,
                   const Zaki::Physics::Date&, const Zaki::String::Directory&) ;
    //.....................................

  //--------------------------------------------------------------
  protected:

    /// virtual 'IClone' method for derived satellites
    virtual Satellite* IClone() const = 0 ;

    /// Exposure time-fraction as a function of energy (GeV)
    virtual double ExpTimeFrac(double energy) const ;

    /// Acceptance (m^2*sr) of the detector as a function of energy (GeV)
    virtual double Acceptance(double energy)  const ;

    /// Numerical scaling so that we can find minos errors with no problems
    /// This will scale the exposure, so it should be multiplied by mu
    /// parameter from the fit
    double NUM_SCALING      = 1e-3 ;
  //--------------------------------------------------------------
  protected:

    std::vector<double> exposure_set ;

    /// Exposure histogram (deep cloned via the 'Clone()' method.)
    TH1F exp_hist ;

    /// Signal shape histogram (deep cloned via the 'Clone()' method.)
    TH1F sig_shape_hist ;

    Data results;

    /// Filed of view (FOV) defined as the half-angle from the normal
    /// to the surface of the detector and in a conical shape.
    double ConeFOV = 90 ; // degrees

    /// Duration of exposure in days
    Zaki::Physics::DateInterval t_duration ;

    /// Actual duration of the data in days
    Zaki::Physics::DateInterval t_duration_actual ;

    /// Inclination of the orbit
    double orbit_inclination      ;

    /// Eccentricity of the orbit
    double eccentricity           ;

    /// Height of the orbit in km
    double orbit_height           ;

    /// Radius of the orbit in km
    double orbit_radius           ;

    /// Period of the orbit in minutes
    double orbit_period           ;

    /// Angular frequency of the orbit (rad/min)
    double orbit_ang_freq         ;

    /// The precession rate of the orbit
    double precession_rate        ;

    /// Right ascension of the ascending node at t_0
    double init_right_asc_node    ;

    /// True anomaly of the orbit at t_0
    double init_true_anomaly      ;

    /// yaw, pitch, roll
    Zaki::Physics::YPR ypr ;

    /// The normal vector to the detector's surface
    ROOT::Math::XYZVector normal  ;

    // Flags for tracking the settings
    bool set_orbit_inclination_flag   = false ;
    bool set_time_duration_flag       = false ;
    // bool set_ecliptic_solar_long_rate_flag = false ;
    bool set_cone_fov_flag        = false ;
    bool set_eccentricity_flag    = false ;
    bool set_orbit_height_flag    = false ;
    bool set_orbit_radius_flag    = false ;
    bool set_orbit_period_flag    = false ;
    bool set_orbit_ang_freq_flag  = false ;
    bool set_right_asc_node_flag  = false ;   
    bool set_init_true_anomaly_flag = false ;
    bool set_YPR_flag             = false ;
    bool set_exposure_eval_flag   = false ;
    bool set_data_flag            = false ;

    bool found_precession_rate_flag = false ;
    bool forced_sun_sync_flag       = false ;
};

//==============================================================
#endif /*DMSS_Satellite_H*/
