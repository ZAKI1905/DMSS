#ifndef DMSS_Analysis_H
#define DMSS_Analysis_H

//  Minuit2
#include <Minuit2/Minuit2Minimizer.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnPrint.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnMinos.h>
#include <Minuit2/MinosError.h>

// Root
#include <TCanvas.h>
#include <TStyle.h>
#include  <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TPaveText.h>

// CONFIND
#include <Confind/ContourFinder.hpp>

// Zaki
#include <Zaki/Math/Math_Core.hpp>

#include "DMSS/Model.hpp"
#include "DMSS/SatBundle.hpp"

//==============================================================
class Analysis : public Prog
{
  //--------------------------------------------------------------
  public:
    
    // Constructor
    Analysis() ;

    // Copy constructor
    Analysis(const Analysis&) ;

    // Assignment operator
    Analysis& operator=(const Analysis&) ;

    // Destructor
    ~Analysis() ;

    // ....................................
    // Setters
    // ....................................

    /// Sets the work directory
    void SetWrkDir(const Zaki::String::Directory&) override;

    /// Sets the model pointer
    void SetModel(std::shared_ptr<Model>) ;

    /// Attaches the satellite (pointer) to the analysis
    void AttachSatellite(std::shared_ptr<Satellite>) ;

    /// Sets the bin periods
    void SetBinPeriod(const size_t& sat_idx, const std::vector<int>&) ;

    /// Sets the energy cut
    void SetECut(double) ;

    /// Sets the confidence level for contour plots
    void SetConfidenceLevel(const double) ;

    /// Sets the plot mode
    void SetPlotMode(const PlotMode&) ;
    // ....................................

    // ....................................
    // Getters
    // ....................................
    /// Gets the bin periods
    std::vector<int> GetBinPeriod(const size_t& sat_idx) const ;

    /// Gets the energy cut
    double GetECut() const ;

    /// Gets the satellite pointer
    std::shared_ptr<Satellite> GetSatellite(const size_t&) ;

    /// Gets the threshold limit fit results
    std::vector<double> GetThreshLimits() const ;

    /// Gets the boosted limit fit results
    double GetBoostLimit() const ;
    // ....................................

    // ....................................
    //   Plotters
    // ....................................

    /// Plots the threshold fit results
    void PlotThresh(const size_t&, const Zaki::String::Directory&) ;

    /// Plots the boosted fit results
    void PlotBoost(const Zaki::String::Directory&) ;
    // ....................................

    //.....................................

    /// Scans the likelihood function vs parameters (Done by Minuit)
    void ScanParMinuit(const size_t& b, const size_t& sat_idx) ;

    /// Scans the likelihood function vs parameters
    std::vector<double> ScanParThresh(const size_t& b_idx, const size_t& t_idx, const size_t& sat_idx) ;

    /// Scans the likelihood function vs parameters
    std::vector<double> ScanParBoost(const double& in_mass, const size_t& sat_idx) ;
    
    /// If the limits are already found and you just want to plot the 
    /// contours in the model's parameter space, set the bypass flag to true
    void PlotParamSpaceThresh(const size_t& sat_idx, const bool bypass=false) ;
    void PlotParamSpaceBoost(const size_t& sat_idx, const bool bypass=false) ;


    /// Performs the threshold limit fit
    /// if no input indices, it will fit all the 
    /// satellite's bins one-by-one.
    void FitThreshold(const size_t& sat_idx, std::vector<size_t>* = nullptr) ;

    /// Fits an individual bin
    void FitBin(Bin, const size_t&) ;

    /// Performs the boosted limit fit
    void FitBoosted(Zaki::Math::Range<double>);

    /// Creates the histograms of fitting each bin (called from 'FitBin')
    double Histgen(const size_t&, 
                   ROOT::Minuit2::MinosError, 
                   ROOT::Minuit2::MinosError, Bin, TH1F);

    // Do analysis
    /// Performs the threshold analysis (Minuit)
    void DoThreshMinuit(const size_t&, std::vector<size_t>* = nullptr, Zaki::Math::Grid2D* = nullptr) ; // Threshold

    /// Performs the threshold analysis
    void DoThresh(const size_t&, std::vector<size_t>* = nullptr, Zaki::Math::Grid2D* = nullptr) ; // Threshold

    /// Function for plotting contours
    /// Not to be confused with the same methods in 'Model'.
    /// Threshold case, (L_dec, G_ann) space
    double ContFuncThresh_LG(double, double)  ;

    /// Threshold case, (M_DM, G_ann) space
    double ContFuncThresh_MG(double, double)  ;


    /// Threshold case, (L_dec, M_DM) space
    double ContFuncThresh_LM(double, double)  ;

    /// Performs the boosted analysis
    void DoBoost(Zaki::Math::Grid2D* = nullptr)  ; // Boosted
    
    /// Function for plotting contours
    /// Not to be confused with the same method in 'Model'.
    /// Boosted case, (L_dec, G_ann) space
    double ContFuncBoost_LG(double, double)   ;
    
    /// Boosted case, (M_DM, G_ann) space
    double ContFuncBoost_MG(double, double)  ;

    /// Boosted case, (L_dec, M_DM) space
    double ContFuncBoost_LM(double, double)  ;

    /// Updates the input Member Function Contour Wrapper (MFCW)
    void UpdateMFCW(CONFIND::MemFuncContWrapper<Analysis, 
                    double (Analysis::*) (double, double)>* ) ;
  //--------------------------------------------------------------
  private:

    bool set_satPtr_flag = false ;
    bool cpy_ctr_called  = false ;
    bool assign_op_called  = false ;

    PlotMode plot_mode = PlotMode::Ldec_Gann ;
    TH1F sig_shape_hist         ;
    std::vector<Bin> focus_bins ;

    /// The confidence level for contour plots:
    double conf_level = 0.95 ;

    std::vector<double> thresh_limits ;
    std::vector<double> thresh_limits_true ;
    std::vector<double> e_bins    ;
    std::vector<int> focus_bin_periods  ;
    std::vector<ROOT::Minuit2::MinosError> boost_fit_results ;
    double boost_fit_val ; 
    double e_cut_val = 50 ; // 50 GeV

    // Pointer member elements
    std::vector<SatBundle> m_SatBundles  ;

    std::shared_ptr<Model> modelPtr = nullptr ;

    void ResetContainers() ;
};

//==============================================================
#endif /*DMSS_Analysis_H*/
