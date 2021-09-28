#include <cassert>

// Root
#include <TF1.h>
#include <TMath.h>
#include <Math/Util.h>
#include <Math/QuantFuncMathMore.h>
#include <gsl/gsl_multiroots.h>

#include "DMSS/LogLikeli.hpp"
#include "DMSS/HybPdf.hpp"

// HybPdf sigbg_hybrid_pdf ;
// HybPdf bg_hybrid_pdf ;

//--------------------------------------------------------------
double ROOT::Minuit2::LogLikeli::operator()(const std::vector<double>& par) const 
{    
  assert(par.size() >= 2) ;

  // Fit parameters
  double mu     = par[0];

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
void ROOT::Minuit2::LogLikeli::AddObsCounts(const Bin& b) 
{
  obs_set.push_back(b.GetTBinObsSet()) ;
  fix_bg = b.GetTBinObsSet()[0].val ;
}

//--------------------------------------------------------------
void ROOT::Minuit2::LogLikeli::AddObsSet
(const std::vector<Zaki::Math::Quantity>& in_obs_set) 
{
  obs_set.push_back(in_obs_set) ;
  fix_bg = in_obs_set[0].val ;
}

//--------------------------------------------------------------
void ROOT::Minuit2::LogLikeli::AddSigShape(const TH1F& sig_shape) 
{
  sig_shape_hist.push_back(sig_shape) ;
}

//--------------------------------------------------------------
// Returns  -2Ln(likelihood_ratio)
double ROOT::Minuit2::LogLikeli::Q(double sig_str, const std::vector<double>& bg_val) const
{
  double val    = 0 ;

  // double penalty = 0 ;

  // Log gets negative when the argument of log approaches '0'
  // double K_penalty = -1e-3 ;

  double sig_bg_rate  = 0 ;

  
  for(size_t i=0 ; i<obs_set.size() ; ++i)
    for (size_t j = 0; j < obs_set[i].size(); j++)
  {

    sig_bg_rate = sig_str*sig_shape_hist[i][ j+1 ] + bg_val[i] ;

    // sigbg_hybrid_pdf.SetPars({sig_bg_rate, obs_set[i][j].val, obs_set[i][j].err}) ;
    // bg_hybrid_pdf.SetPars({bg_val[i], obs_set[i][j].val, obs_set[i][j].err}) ;
    
    // sig_bg = sigbg_hybrid_pdf.Integrate() ;
    // bg     = bg_hybrid_pdf.Integrate() ;

    // if (obs_set[i][j].val > 50)
    // {
    //   sig_bg = TMath::Gaus(obs_set[i][j].val, sig_bg_rate, sqrt(sig_bg_rate)) ;
    //   // bg =  TMath::Gaus(obs_set[i][j].val, bg_val[i], sqrt(bg_val[i])) ;
    // }
    // else
    // {
    //   sig_bg = TMath::Poisson(obs_set[i][j].val, sig_bg_rate) ;
    //   // bg = TMath::Poisson(obs_set[i][j].val, bg_val[i]) ;
    // }

    /// safe evaluation of log(x) with a protections against negative or zero argument 
    val += obs_set[i][j].val*ROOT::Math::Util::EvalLog(sig_bg_rate) - sig_bg_rate ;


    // penalty += K_penalty*ROOT::Math::Util::EvalLog(bg_val[i] - 0.1) ;
    // penalty += K_penalty*ROOT::Math::Util::EvalLog(bg_val[i] - 1000) ;
  }
  // std::cout << " -> " << val << "\n" ;
  return -2*val ;
}

//--------------------------------------------------------------
/// -2LogL = 5.99146 Equation
double ROOT::Minuit2::LogLikeli::MuEquation(double x)
{
  return this->operator()({x, var_bg}) 
          - this->operator()({1e-7, fix_bg}) 
          - ROOT::MathMore::chisquared_quantile(0.95, 2) ;
}

//--------------------------------------------------------------
/// F'(x)
double ROOT::Minuit2::LogLikeli::MuEquationDer(double x)
{
  double val  = 0 ;
  double tmp  = 0 ;

  for(size_t i=0 ; i<obs_set.size() ; ++i)
    for (size_t j = 0; j < obs_set[i].size(); j++)
  {

    tmp = sig_shape_hist[i][ j+1 ] / (x * sig_shape_hist[i][ j+1 ] + var_bg) ;
    // tmp = sig_shape_hist[i][ j+1 ] / (x * sig_shape_hist[i][ j+1 ] + bg_val[i]) ;

    /// safe evaluation of log(x) with a protections against negative or zero argument 
    val += obs_set[i][j].val * tmp -  sig_shape_hist[i][ j+1 ] ;
  }

  return -2*val ;
}
//--------------------------------------------------------------
/// { F(x), F'(x) }
void ROOT::Minuit2::LogLikeli::MuEquationFdf(double x, double* in_f, double* in_df) 
{
  *in_f   = this->operator()({x, var_bg}) 
          - this->operator()({1e-7, fix_bg}) 
          - ROOT::MathMore::chisquared_quantile(0.95, 2) ;


  
  double val  = 0 ;
  double tmp  = 0 ;

  for(size_t i=0 ; i<obs_set.size() ; ++i)
    for (size_t j = 0; j < obs_set[i].size(); j++)
  {

    tmp = sig_shape_hist[i][ j+1 ] / (x * sig_shape_hist[i][ j+1 ] + var_bg) ;
    // tmp = sig_shape_hist[i][ j+1 ] / (x * sig_shape_hist[i][ j+1 ] + bg_val[i]) ;

    /// safe evaluation of log(x) with a protections against negative or zero argument 
    val += obs_set[i][j].val * tmp -  sig_shape_hist[i][ j+1 ] ;
  }        

  *in_df  = -2*val ;
}

//--------------------------------------------------------------
void ROOT::Minuit2::LogLikeli::SetVarBg(double in_var_bg)
{
  var_bg = in_var_bg ;
}

//--------------------------------------------------------------
double ROOT::Minuit2::LogLikeli::ThreshCont(double mu, double bg)
{
  return this->operator()({mu, bg}) - this->operator()({1e-7, fix_bg}) ;
}
//--------------------------------------------------------------
void ROOT::Minuit2::LogLikeli::SetNDBgSet(const size_t ignore_idx, const std::vector<double>& in_nd_bg_set) 
{
  nd_bg_set = in_nd_bg_set ;
  nd_ignore_idx = ignore_idx ;
}
//--------------------------------------------------------------
void ROOT::Minuit2::LogLikeli::SetNDBestFit(const std::vector<double>& in_bestfit) 
{
  nd_best_fit = in_bestfit ;
}
//--------------------------------------------------------------
double ROOT::Minuit2::LogLikeli::BoostCont(double mu, double bg)
{
  std::vector<double> par = {mu} ;
  for (size_t i = 0; i < nd_bg_set.size() ; i++)
  {
    if( i != nd_ignore_idx)
      par.push_back(nd_bg_set[i]) ;
    else
    {
      par.push_back(bg) ;
    }
  }
  
  return this->operator()(par) - this->operator()(nd_best_fit) ;
}
//--------------------------------------------------------------

/// -2LogL = 5.99146 Equation
// double ROOT::Minuit2::LogLikeli::MuEquationND(const gsl_vector* x, gsl_vector* in_f)
// {
//   const double x0 = gsl_vector_get(x,0);
//   const double x1 = gsl_vector_get(x,1);

//   const double y0 = pars.a * x0  - 4 ;
//   const double y1 = pars.b * x0 + x1 * x0 ;

  
//   return this->operator()({x, fix_bg}) 
//           - this->operator()({1e-7, fix_bg}) 
//           - ROOT::MathMore::chisquared_quantile(0.95, 2) ;
// }

//--------------------------------------------------------------


