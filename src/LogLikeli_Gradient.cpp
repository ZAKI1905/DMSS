#include <cassert>

// Root
#include <TF1.h>
#include <TMath.h>
#include <Math/Util.h>

// #include <Zaki/Math/GSLFuncWrapper.hpp>
#include <Zaki/Math/Integrator.hpp>

#include "DMSS/LogLikeli_Gradient.hpp"
#include "DMSS/HybPdf.hpp"

HybPdf sigbg_hybrid_pdf ;
HybPdf bg_hybrid_pdf ;

//--------------------------------------------------------------
namespace Zaki::Math 
{
//--------------------------------------------------------------
class LogLikeliGradIntegral : public Integrator
{

  protected:
    double rate  ;
    Quantity obs ;

  public:
    LogLikeliGradIntegral() : Integrator(1000, 1e-10, 1e-10) {}
    virtual ~LogLikeliGradIntegral() {}
    virtual Integrator* GetIntegratorPtr() = 0 ;
    virtual double Integrand(double) const = 0 ;

    void SetParameters(const double& in_rate, const Quantity&) ;

};

  // Integrator* LogLikeliGradIntegral::GetIntegratorPtr()
  // {
  //   return this;
  // }

  // double LogLikeliGradIntegral::Integrand(double x) const
  // {
  //   return TMath::Poisson(x-1, rate)*TMath::Gaus(x, obs.val, obs.err) ;
  // }

void LogLikeliGradIntegral::SetParameters(const double& in_rate, const Quantity& in_obs)
{
  rate = in_rate ;
  obs = in_obs ;
}
//--------------------------------------------------------------
class NumeratorIntegral : public LogLikeliGradIntegral
{

  public:
    NumeratorIntegral() {}
    ~NumeratorIntegral() {}
    
    virtual NumeratorIntegral* GetIntegratorPtr() override
    {
      return this;
    }

    virtual double Integrand(double x) const override
    {
      return TMath::Poisson(x-1, rate)*TMath::Gaus(x, obs.val, obs.err, true) ;
    }
};
//--------------------------------------------------------------
class DenomIntegral : public LogLikeliGradIntegral
{

  public:
    DenomIntegral() {}
    ~DenomIntegral() {}
    
    virtual DenomIntegral* GetIntegratorPtr() override
    {
      return this;
    }

    virtual double Integrand(double x) const override
    {
      return TMath::Poisson(x, rate)*TMath::Gaus(x, obs.val, obs.err, true) ;
    }
};
//--------------------------------------------------------------
}

Zaki::Math::NumeratorIntegral num_int ;
Zaki::Math::DenomIntegral denom_int ;

//--------------------------------------------------------------
double ROOT::Minuit2::LogLikeliGradient::operator()(const std::vector<double>& par) const 
{    
  assert(par.size() >= 2) ;

  // Fit parameters
  double mu     = par[0];
  // double bg_val = par[1];

  std::vector<double>::const_iterator first = par.begin() + 1;
  std::vector<double>::const_iterator last = par.end();
  std::vector<double> bg_val(first, last);

  // std::cout << "par[0]: "<< par[0] << " , ";
  // std::cout << "par[1]: "<< par[1] << " , ";
  // std::cout << "par[2]: " << par[2] << "\n";


  double Q_val = 0 ; 

  Q_val +=  Q(mu, bg_val);

  return Q_val;
}

//--------------------------------------------------------------
void ROOT::Minuit2::LogLikeliGradient::AddObsCounts(const Bin& b) 
{
  obs_set.push_back(b.GetTBinObsSet()) ;
}

//--------------------------------------------------------------
void ROOT::Minuit2::LogLikeliGradient::AddSigShape(const TH1F& sig_shape) 
{
  sig_shape_hist.push_back(sig_shape) ;
}

//--------------------------------------------------------------
// Returns  -2Ln(likelihood_ratio)
double ROOT::Minuit2::LogLikeliGradient::Q(double sig_str, const std::vector<double>& bg_val) const
{
  double val    = 0 ;

  double penalty = 0 ;

  // Log gets negative (-inf) when the argument of log approaches '0'
  double K_penalty = -1e-2 ;

  // HybPdf func;

  double sig_bg, bg, sig_bg_rate  = 0 ;

  
  for(size_t i=0 ; i<obs_set.size() ; ++i)
    for (size_t j = 0; j < obs_set[i].size(); j++)
  {

    // double pars[] = {sig_str*sig_shape_hist[i][ j+1 ] + bg_val[i], 
    //                   obs_set[i][j].val, obs_set[i][j].err};

    // func.SetPars(pars) ;

    // if(bg_val[i] < 1 || bg_val[i]> 500)
    // {
    //   val -= 99999999 ;
    //   continue ;
    // }


    sig_bg_rate = sig_str*sig_shape_hist[i][ j+1 ] + bg_val[i] ;

    sigbg_hybrid_pdf.SetPars({sig_bg_rate, obs_set[i][j].val, obs_set[i][j].err}) ;
    bg_hybrid_pdf.SetPars({bg_val[i], obs_set[i][j].val, obs_set[i][j].err}) ;
    
    sig_bg = sigbg_hybrid_pdf.Integrate() ;
    bg     = bg_hybrid_pdf.Integrate() ;

    // if (obs_set[i][j].val > 50)
    // {
    //   sig_bg = TMath::Gaus(obs_set[i][j].val, rate, sqrt(rate)) ;
    //   bg =  TMath::Gaus(obs_set[i][j].val, bg_val[i], sqrt(bg_val[i])) ;
    //   // std::cout<< "====>> obs_set["<<i<<"]["<<j<<".val  "<< obs_set[i][j].val
    //   // << ", bg_val[i]: " << bg_val[i] << ", sig_shape_hist[i][j+1]: "
    //   // << sig_shape_hist[i][ j+1 ] 
    //   // << ",  rate: " << rate <<"\n" ;

    // }
    // else
    // {
      // sig_bg = TMath::Poisson(obs_set[i][j].val, sig_bg_rate) ;
      // bg = TMath::Poisson(obs_set[i][j].val, bg_val[i]) ;
    // }

    // if(bg_val[i] <= 1)
    // {
    
    //   std::cout <<" bg= "<< bg << ", bg_val["<<i<<"]= " << bg_val[i] << ", sig_bg_rate= "<< sig_bg_rate << ", sig_bg=" << sig_bg <<"\n" ;
    // }
    /// safe evaluation of log(x) with a protections against negative or zero argument 
    val += ROOT::Math::Util::EvalLog(sig_bg / bg) ;

    // penalty += K_penalty*ROOT::Math::Util::EvalLog(bg_val[i] - 0.499) ;
    // penalty += K_penalty*ROOT::Math::Util::EvalLog(200.1 - bg_val[i]) ;
  }
  
  return -2*val + penalty ;
}

//--------------------------------------------------------------
std::vector<double> ROOT::Minuit2::LogLikeliGradient::Gradient(const std::vector<double>& in_pars) const 
{
  // The gradient vector 
  // 1-st component is the signal strength mu,
  // the rest are background pars.
  // Note that: obs_set.size() = in_pars.size() by construction
  std::vector<double> del(obs_set.size(), 0) ;
  
  // Loop over the energy bins
  for(size_t i=0 ; i < obs_set.size() ; ++i)
  {
    // The assumption is that the observed events are the same in all 
    // time-bins of the same energy bin, i.e. all the same i's but 
    // different j's which is why we perform half the integrals (I_0)
    // in the first loop to save time.
    double I_0 = GetIntRatio(in_pars[i+1], obs_set[i][0]) ;

    // Loop over the time-bins
    for (size_t j = 0; j < obs_set[i].size(); j++)
    {

      // Note that: obs_set[i][j] = obs_set[i][0] 
      double I_1 = GetIntRatio(in_pars[0]*sig_shape_hist[i][j+1] + in_pars[i+1],
                               obs_set[i][j]) ;

      // This is the mu component derivate:
      del[0] += -2*sig_shape_hist[i][ j+1 ]*(I_1 - 1.0);

      // This is the bg_i component derivate:
      del[i] += -2*(I_1 - I_0);
    }
  }

  return del ;

}
//--------------------------------------------------------------
double ROOT::Minuit2::LogLikeliGradient::GetIntRatio(const double& rate, 
                                  const Zaki::Math::Quantity& in_obs) const
{
  num_int.SetParameters(rate, in_obs) ;
  num_int.Integrate(in_obs.val - in_obs.err,
                in_obs.val + in_obs.err ) ;

  denom_int.SetParameters(rate, in_obs) ;
  denom_int.Integrate(in_obs.val - in_obs.err,
                      in_obs.val + in_obs.err ) ;

  return (num_int.GetIntegrationResult().val / denom_int.GetIntegrationResult().val) ;
}
//--------------------------------------------------------------

//--------------------------------------------------------------

