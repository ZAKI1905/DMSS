#ifndef DMSS_Fitter_H
#define DMSS_Fitter_H

#include "DMSS/Satellite.hpp"

//==============================================================
class Fitter : public Prog
{
  public:
    
    /// Default constructor
    Fitter() ;

    /// Copy constructor
    Fitter(const Fitter&) = delete;

    /// Destructor
    ~Fitter() ;

    /// Perform the fitting
    std::vector<double> Fit() ;
    double Fit2D() const ;


    void AddObsCounts(const Bin&) ;
    void AddSigShape(const TH1F&) ;

  private:
    // double theErrorDef = 1 ;
    double FitIdx(const size_t&) const ;
    std::vector<std::vector<Zaki::Math::Quantity> > obs_set ;
    std::vector<TH1F> sig_shape_hist ;
    std::vector<double> fit_results ;

};



//==============================================================
#endif /*DMSS_Fitter_H*/
