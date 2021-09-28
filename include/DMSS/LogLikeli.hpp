#ifndef DMSS_LogLikeli_H
#define DMSS_LogLikeli_H

#include <Minuit2/FCNBase.h>

#include "DMSS/Bin.hpp"

//==============================================================
namespace ROOT {

namespace Minuit2 {


class LogLikeli : public FCNBase, public Prog
{

  public:

    ///  Constructor
    LogLikeli() : Prog("LogLikeli", true) {}
    
    /// Copy constructor
    LogLikeli(const LogLikeli& other)
        : Prog("LogLikeli", true),
            theErrorDef(other.theErrorDef),
            obs_set(other.obs_set), 
            fix_bg(other.fix_bg), 
            var_bg(other.var_bg),
            nd_bg_set(other.nd_bg_set),
            nd_ignore_idx(other.nd_ignore_idx),
            nd_best_fit(other.nd_best_fit)
            {
                for (auto h: other.sig_shape_hist) {
                    sig_shape_hist.emplace_back( *(static_cast<TH1F*>(h.Clone())) ) ;
                }
            }
    
    virtual double operator()(const std::vector<double>& par) const override ;
    
    ~LogLikeli() {}

    double Up() const override {return theErrorDef;}

    void SetErrorDef (double def) override
    {  
      theErrorDef = def; 
      char tmp[100] ;
      sprintf(tmp, "%.3f", def) ;
      Z_LOG_INFO("The error definition is set to: " + std::string(tmp) + ".") ;
    }
    
    void AddObsCounts(const Bin&) ;
    void AddObsSet (const std::vector<Zaki::Math::Quantity>& in_obs_set) ;
    
    void AddSigShape(const TH1F&) ;
    void SetNDBgSet(const size_t ignore_idx, const std::vector<double>&) ;
    void SetNDBestFit(const std::vector<double>& in_bestfit) ;

    double Q(double, const std::vector<double>&) const;

    /// -2LogL = 5.99146 Equation
    double MuEquation(double x) ;
    double MuEquationDer(double x) ;
    void MuEquationFdf(double x, double* in_f, double* in_df)  ;
    void SetVarBg(double in_var_bg) ;
    double fix_bg, var_bg ;
    double ThreshCont(double mu, double bg) ;
    double BoostCont(double mu, double bg) ;

    
  private:
    double theErrorDef = 1 ;
    std::vector<std::vector<Zaki::Math::Quantity> > obs_set ;
    std::vector<TH1F> sig_shape_hist ;

    std::vector<double> nd_bg_set ;
    size_t nd_ignore_idx  ;
    std::vector<double> nd_best_fit ;
};

  }  // namespace Minuit2

}  // namespace ROOT
//==============================================================
#endif /*DMSS_LogLikeli_H*/
