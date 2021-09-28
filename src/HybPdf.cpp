#include <cassert>

// #include <gsl/gsl_integration.h>

// Root
#include <TMath.h>

#include <Zaki/Math/GSLFuncWrapper.hpp>

// Local Headers
#include "DMSS/HybPdf.hpp"

//--------------------------------------------------------------
HybPdf::HybPdf() 
{   
  GSL_workspace = gsl_integration_workspace_alloc(2000);
}

//--------------------------------------------------------------
HybPdf::~HybPdf() 
{   
  gsl_integration_workspace_free(GSL_workspace); 
}

//--------------------------------------------------------------
void HybPdf::SetPars(const std::vector<double>& par) 
{
  assert(par.size() >= 3) ;

  rate    = par[0] ;
  obs     = par[1] ;
  obs_err = par[2] ;
}

//--------------------------------------------------------------
double HybPdf::IntegrandP(double x)
{
  return TMath::Poisson(x, rate)*TMath::Gaus(x, obs, obs_err, true) ;
}

//--------------------------------------------------------------
double HybPdf::IntegrandG(double x)
{
  return TMath::Gaus(x, rate, sqrt(rate))*TMath::Gaus(x, obs, obs_err, true) ;
}

//--------------------------------------------------------------
double HybPdf::NormalFactor(double x)
{
  return TMath::Gaus(x, obs, obs_err, true) ;
}

//--------------------------------------------------------------
double HybPdf::Integrate()
{
  // ROOT::Math::WrappedMemFunction<HybPdf, double ( HybPdf::* ) (double)>
  //  f(*this, &HybPdf::Integrand);

  // ROOT::Math::WrappedMemFunction<HybPdf, double ( HybPdf::* ) (double)>
  //  g(*this, &HybPdf::NormalFactor);

  // ROOT::Math::Integrator f_int(f);
  // ROOT::Math::Integrator g_int(g);

  double num_err, norm_err;

  // DarkPhoton local_model ;
  // local_model.Init() ;

  Zaki::Math::GSLFuncWrapper<HybPdf, double (HybPdf::*)(double)> f ;

  if(rate > 100) 
    f.SetMemberFunc(this, &HybPdf::IntegrandG); 
  else
    f.SetMemberFunc(this, &HybPdf::IntegrandP); 

     
  Zaki::Math::GSLFuncWrapper<HybPdf, double (HybPdf::*)(double)> 
    g(this, &HybPdf::NormalFactor);     

  gsl_function f_GSL = static_cast<gsl_function> (f) ; 
  gsl_function g_GSL = static_cast<gsl_function> (g) ; 


  double obs_min = ( obs-obs_err > 0 ) ? obs-obs_err : 0 ;

  double num  = 0 ; // f_int.Integral(obs_min, obs+2*obs_err);
  double norm = 1 ; // g_int.Integral(obs_min, obs+2*obs_err);

  gsl_integration_qag(&f_GSL, obs_min, obs+obs_err, 1e-13, 1e-13, 1000, 1,
                      GSL_workspace, &num, &num_err);
  gsl_integration_qag(&g_GSL, obs_min, obs+obs_err, 1e-13, 1e-13, 1000, 1, 
                      GSL_workspace, &norm, &norm_err);

  return num/norm ;
}

//--------------------------------------------------------------
