/* 

  Assuming a constant background, and a randomly 
   (const) generated data, and an oscillating signal,
   we find the 95% limit on the signal strength.
   We vary the time-binning and check how the limit changes.

  - Last updated by Zaki Nov 23, 2020

*/

// #include <Zaki/File/VecSaver.hpp>
#include <Zaki/Util/ObjObserver.hpp>
#include <Zaki/File/VecSaver.hpp>

// Local headers
#include "DMSS/Bin.hpp"
#include "DMSS/LogLikeli.hpp"
// #include "DMSS/GenericModel.hpp"
// #include "DMSS/DarkPhoton.hpp"
// #include "DMSS/Prog.hpp"

#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TRandom3.h>
// #include <TRandom.h>

#include <Math/QuantFuncMathMore.h>

//  Minuit2
#include <Minuit2/Minuit2Minimizer.h>
// #include "Minuit2/MnUserParameters.h"
// #include "Minuit2/VariableMetricMinimizer.h"
#include <Minuit2/FunctionMinimum.h>
// #include "Minuit2/MnUserParameterState.h"
#include <Minuit2/MnPrint.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnMinos.h>
// #include "Minuit2/MnContours.h"
// #include "Minuit2/MnPlot.h"
#include <Minuit2/MinosError.h>
// #include "Minuit2/ContoursError.h"


//************//
//    MAIN
//************//
int main() {
  
  using namespace Zaki::Util ;


  Zaki::String::Directory dir(__FILE__) ;


  ObjManager::SetFile(dir.ParentDir() +"/results/Time_Binning/Objs_time_binning.txt") ;
  LogManager::SetLogLevels(LogLevel::Info, LogLevel::Verbose) ;
  LogManager::SetLogFile(dir.ParentDir() +"/results/Time_Binning/log_time_binning.txt") ;
  // LogManager::SetBlackWhite(true) ;


  size_t n_bins = 100 ;
  double sig_period = 35 ;

for (size_t j = 0; j < 100; j++)
{
  gRandom = new TRandom3(0);
  //......................................................
  // Data
  TH1F h_obs("Events", "Events", n_bins, 0, 100) ;

  TF1 *mygaus = new TF1("flat","1", 0, 100);
  h_obs.FillRandom("flat", 1000) ;
  // hh.Scale( 1.0/hh.GetNormFactor() ) ;
  
  // std::vector<Zaki::Math::Quantity> two_bin_obs_set = { {0, 0}, {0, 0} };

  // for (size_t i = 0; i < n_bins; i += 2)
  // {
  //   two_bin_obs_set[0].val += obs_set[i].val ;
  //   two_bin_obs_set[1].val += obs_set[i+1].val ;
  // }
  //......................................................

  //......................................................
  // Signal
  TH1F hsig("signal", "sig", n_bins, 0, 100) ;

  TF1 *sinsig = new TF1("sinsig","sin(2*pi*x/[0]) + 1", 0, 100);

  sinsig->SetParameter(0, sig_period);

  for (size_t i = 0; i < 100; i++)
  {
    hsig.AddBinContent(i+1, sinsig->Integral(i, i+1)) ;
  }
  // hsig.FillRandom("sinsig", 1000) ;
  hsig.Scale(10) ;



  // TH1F two_bin_sig("signal_2b", "sig2b", 2, 0, 100) ;
  // std::pair<double, double> tmp_sig_2_bins = {0, 0} ;
  // for (size_t i = 0; i < n_bins; i += 2)
  // {
  //   tmp_sig_2_bins.first += hsig.GetBinContent(i+1) ;
  //   tmp_sig_2_bins.second += hsig.GetBinContent(i+2) ;
  // }
  // two_bin_sig.SetBinContent(1, tmp_sig_2_bins.first) ;
  // two_bin_sig.SetBinContent(2, tmp_sig_2_bins.second) ;
  //......................................................

  //......................................................
  // Plotting
  TCanvas c1("c1", "Sin-Flat", 1000, 800) ;

  hsig.SetLineColor(kRed) ;
  hsig.Draw() ; 
  
  h_obs.Draw("same") ;

  c1.SaveAs((dir.ParentDir() +"/results/Time_Binning/T_"+ std::to_string((int)sig_period) + "/Sin_Flat_"
              +std::to_string(n_bins)+".pdf").Str().c_str()) ;
  //......................................................
  std::vector<std::vector<double> > fit_results; fit_results.reserve(50) ;
  std::vector<double> bin_period_set = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
                                        13, 14, 15, 16, 17, 18, 19, 21, 22,
                                        23, 24, 25, 30, 40, 50, 100} ;

  for (size_t i = 0; i < bin_period_set.size(); i++)
  {
    
    double bin_period = bin_period_set[i] ;

    TH1F sig_shape_hist = *(TH1F*)(hsig.Clone(("sig_shape_"+std::to_string(bin_period)).c_str()));
    sig_shape_hist.Rebin(bin_period);

    double scale_factor = 1.0 / sig_shape_hist.Integral();
    sig_shape_hist.Scale(100*scale_factor) ;


    TH1F h_obs_rebin = *(TH1F*)(h_obs.Clone(("obs_"+std::to_string(bin_period)).c_str()));
    h_obs_rebin.Rebin(bin_period) ;

    std::vector<Zaki::Math::Quantity> obs_set ; 
    for (size_t i = 0; i < (1.0*n_bins) / bin_period; i++)
    {
      obs_set.push_back({h_obs_rebin.GetBinContent(i+1), h_obs_rebin.GetBinError(i+1)}) ;
    }

    // Bin b({100, 200}, {150, 5}, {},) ;  
    // size_t chops = hsig.GetSize() -2  ;
    // b.divide(chops);

    ROOT::Minuit2::LogLikeli  fcn;

    fcn.AddObsSet(obs_set) ;
    fcn.AddSigShape(sig_shape_hist)  ;


    // fcn.AddObsSet(two_bin_obs_set) ;
    // fcn.AddSigShape(two_bin_sig) ;
    //......................................................
    //  Minimize
    //......................................................

    //(para_name, para_start_value, para_uncertainties)
    ROOT::Minuit2::MnUserParameters upar ;
    upar.Add("mu", 1e-1, 1e-1);
    upar.Add("bg", bin_period*(1000.0/n_bins), 2);

    // set limits
    upar.SetLowerLimit("mu", 0);
    upar.SetLowerLimit("bg", 0);

    // Creat MIGRAD minimizer
    int strategy_counter = 2 ;
    ROOT::Minuit2::MnMigrad migrad(fcn, upar, strategy_counter);

    // Minimize
    ROOT::Minuit2::FunctionMinimum min = migrad();

    // output
    std::cout<<"minimum: "<<min<<std::endl;

    //......................................................
    // ****************
    //  Error analysis
    // ****************
    // MINOS Error analysis
    ROOT::Minuit2::MnMinos Minos(fcn, min);

    // MINOS errors
    // -2Log & chi2 have the same normalization:
    // Zaki::Math::NDimContLevel cont_level((int)upar.Params().size(), conf_level) ;
    // fcn.SetErrorDef(cont_level.Solve().Up);
    fcn.SetErrorDef(ROOT::MathMore::chisquared_quantile(0.95,
                    (int)upar.Params().size())) ;

    ROOT::Minuit2::MinosError e0 = Minos.Minos(0); // get parameter 0's error
    ROOT::Minuit2::MinosError e1 = Minos.Minos(1); // get parameter 1's error

    // output
    Z_LOG_INFO("Minos errors: ") ;
    std::cout<<e0<<"\n"<<e1<<"\n";

    char tmp_char[150] ;
    sprintf(tmp_char, "\n -Lower Limit (95%%): %.2e\n -Upper Limit (95%%): %.2e",
            e0.Min() + e0.Lower(), e0.Min() + e0.Upper()) ;
    Z_LOG_INFO(tmp_char) ;

    if(!e0.IsValid() || !e1.IsValid() ) Z_LOG_ERROR(" ---> Minos error is invalid!") ;
    //......................................................

    //......................................................
    // Plotting
    TCanvas c2("c2", "Sin-Flat Fit", 1000, 800) ;

    sig_shape_hist.Scale(e0.Min() + e0.Upper()) ;

    TH1F sig_bg_hist("signal-bg", "sig-bg", (size_t) ((1.0*n_bins) / bin_period), 0, 100) ;
    for (size_t i = 0; i < (size_t) ((1.0*n_bins) / bin_period); i++)
    {
      sig_bg_hist.SetBinContent(i+1, e1.Min() + sig_shape_hist.GetBinContent(i+1)) ;
      // sig_bg_hist.SetBinError(i+1, e1.Upper()) ;
    }

    h_obs_rebin.Draw() ;
    
    h_obs_rebin.SetMaximum(1.05*std::max(h_obs_rebin.GetMaximum(), sig_bg_hist.GetMaximum()) ) ;
    h_obs_rebin.SetMinimum(0.95*std::min(h_obs_rebin.GetMinimum(), sig_bg_hist.GetMinimum()) ) ;
    
    sig_bg_hist.SetLineColor(kRed) ;
    sig_bg_hist.Draw("same") ;


    
    c2.SaveAs((dir.ParentDir() +"/results/Time_Binning/T_"
              + std::to_string((int)sig_period) + "/Sin_Flat_Fit_"
              + std::to_string((int)bin_period)+".pdf").Str().c_str()) ;
    //......................................................

    fit_results.push_back({bin_period, e0.Min() + e0.Upper()}) ;
  }

  Zaki::File::VecSaver saver(dir.ParentDir() +"/results/Time_Binning/T_"
                              + std::to_string((int)sig_period) 
                              + "/Sin_Flat_Fit_Results_"+std::to_string(j)+".tsv") ;

  saver.Export2D(fit_results, "\t") ;

}
  return 0;
}
