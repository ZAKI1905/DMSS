#ifndef DMSS_DarkPhoton_H
#define DMSS_DarkPhoton_H

// Zaki
#include <Zaki/Physics/Constants.hpp>

#include "DMSS/Data.hpp"
#include "DMSS/Model.hpp"

//==============================================================

//==============================================================
// Dark photon model parameters
struct DPhotPar
{

  double mDM; // Dark matter mass
  double mDP; // Dark photon mass
  double eps; // Kinetic mixing of A-A'
  double gX ; // DM-Dark photon gauge interaction coupling

  double alphaX() const { return pow(gX,2)/4/M_PI ; }; // Structure constant
};

//==============================================================
class DarkPhoton : public Model
{

  //--------------------------------------------------------------
  public:

    DarkPhoton()  ;
    ~DarkPhoton() ;

    //............................................
    // Setters
    //............................................
    void SetModelPars(DPhotPar&);
    void SetDMMass(const double) override ;
    void SetDMMass(const std::vector<double>&) override ;
    void SetElemSet();

    void FixRelic() ;

    void Init() override ;

    // Importing the dark photon branching ratio
    void ImportDarkPhotonBr(const Zaki::String::Directory&) ;

    //............................................

    //............................................
    // Getters
    //............................................
    // std::string GetLabel() ;
    DPhotPar GetModelPars() const;

    // Velocity of thermal DM in the sun
    double GetSunDMVel()    const ;

    // Sommerfeld Enhancement Factor
    double GetSommerfeld() const ;

    // Sommerfeld (with resonances)
    // add later

    // DM annihilation cross section 
    double GetAnnXSec() const ;
    double GetCAnn() const ;
    double GetEqTau() const ;
    double GetAnnRate() const ;

    // Elastic scattering for a specific element
    double GetElasticSig(Zaki::Physics::Element) const ;

    double GetCapRate(Zaki::Physics::Element) const ;
    double GetCapRateTot() const ;

    // Dark photon decay lengths
    double GetDecLenT() const ; // Threshold Case
    double GetDecLenB() const ; // Boosted Case

    // Dark photon branching ratio
    double GetDPhotonBr() const ;
    double DarkPhotonBrWrapper(double*, double*) const ;

    // Dark photon decay probability (threshold)
    double GetDecayProbT() const ;
    double GetDecayProbB() const ;
    //............................................
    

    //............................................
    // Plotters
    //............................................
    // Plot 
    // void Plot() const;

    // Overriding the base class (Model) method
    double ContFuncThresh_LG(double, double) override ;
    double ContFuncThresh_MG(double, double) override ;
    double ContFuncThresh_LM(double, double) override ;

    // Boosted case
    double ContFuncBoost_LG(double, double) override ;
    double ContFuncBoost_MG(double, double) override ;
    double ContFuncBoost_LM(double, double) override ;

    // Plot dark photon branching ratio
    void DPhotonBrPlot(double, double, const std::string&) const;
    double SommIntegrand(double) ;

    
  //--------------------------------------------------------------
  private:

    virtual DarkPhoton* IClone() const override ;
    
    double sommerfeld_factor = 1;

    std::vector<double>  br_vals ;
    std::vector<double>  m_vals  ;

    DPhotPar dark_model ;
    std::vector <Zaki::Physics::Element> element_set ;
    std::vector<std::pair<double, double> > DP_mass_br ;

    // Flags for tracking the settings
    bool set_model_par_flag = false ;
    bool fix_relic_flag = false ;
    bool set_elem_set_flag = false ;
    bool set_mass_br_flag = false ;
    bool exact_somm_flag  = true ;
    bool found_exact_somm_flag = false ;
    
    // Approximate Sommerfeld Enhancement formula (no resonances)
    double GetApproxSomm()  const ;

    // Exact Sommerfeld Enhancement
    double GetExactSomm()  const ;

    // Finds Sommerfeld Enhancement Factor
    void EvalExactSomm() ;

};

//==============================================================
#endif /*DMSS_DarkPhoton_H*/
