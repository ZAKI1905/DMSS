#include <cassert>

#include "DMSS/ChiSqrd.hpp"


//--------------------------------------------------------------
double ROOT::Minuit2::ChiSqrd::operator()(const std::vector<double>& par) const
{
  assert(par.size() >= 2);

  // Fit parameters
  double mu     = par[0];
  double bg_val = par[1];

  double predict_val = 0 ;

  double chi2_val = 0 ;

  for(size_t i=0 ; i<obs_set.size() ; ++i) 
  {
    for (size_t j = 0; j < obs_set[i].size(); j++)
    {
      predict_val = mu*sig_shape_hist[i][ j+1 ] + bg_val;
      chi2_val +=  pow(obs_set[i][j].val - predict_val, 2) / pow(obs_set[i][j].err, 2);
    }
  }

  return chi2_val;
}

//--------------------------------------------------------------
void ROOT::Minuit2::ChiSqrd::AddObsCounts(const Bin& b) 
{
  obs_set.push_back(b.GetTBinObsSet()) ;
}

//--------------------------------------------------------------
void ROOT::Minuit2::ChiSqrd::AddSigShape(const TH1F& sig_shape) 
{
  sig_shape_hist.push_back(sig_shape) ;
}

//--------------------------------------------------------------