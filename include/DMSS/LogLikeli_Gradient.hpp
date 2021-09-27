#ifndef DMSS_LogLikeli_Gradient_H
#define DMSS_LogLikeli_Gradient_H

#include <Minuit2/FCNGradientBase.h>

#include "DMSS/Bin.hpp"

//==============================================================
namespace ROOT {

namespace Minuit2 {


class LogLikeliGradient : public FCNGradientBase, public Prog
{

  public:

    ///  Constructor
    LogLikeliGradient() : Prog("LogLikeliGradient", true) {}
    
    /// Copy constructor
    LogLikeliGradient(const LogLikeliGradient& other)
        : Prog("LogLikeliGradient", true),
            theErrorDef(other.theErrorDef),
            obs_set(other.obs_set)
            {
                for (auto b: sig_shape_hist) {
                    b = *(static_cast<TH1F*>(b.Clone())) ;
                }
            }
    
    virtual double operator()(const std::vector<double>& par) const override ;
    
    ~LogLikeliGradient() {}

    double Up() const override {return theErrorDef;}

    void SetErrorDef (double def) override
    {  
      theErrorDef = def; 
      char tmp[100] ;
      sprintf(tmp, "%.3f", def) ;
      Z_LOG_INFO("The error definition is set to: " + std::string(tmp) + ".") ;
    }

   virtual bool CheckGradient() const override {return false;}

    virtual std::vector<double> Gradient(const std::vector<double>&) const override ;
    double GetIntRatio(const double& rate, const Zaki::Math::Quantity& in_obs) const ;

    void AddObsCounts(const Bin&) ;
    void AddSigShape(const TH1F&) ;
    double Q(double, const std::vector<double>&) const;

  private:
    double theErrorDef = 1 ;
    std::vector<std::vector<Zaki::Math::Quantity> > obs_set ;
    std::vector<TH1F> sig_shape_hist ;
};

  }  // namespace Minuit2

}  // namespace ROOT
//==============================================================
#endif /*DMSS_LogLikeli_Gradient_H*/
