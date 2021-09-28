/*
  Fitter class

*/

#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooHistPdf.h>
#include <RooAddPdf.h>
#include <RooPlot.h>
#include "RooFitResult.h"
#include "RooMinimizer.h"
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>

#include "DMSS/Fitter.hpp"
//==============================================================
//--------------------------------------------------------------
// Constructor
Fitter::Fitter() : Prog("Fitter")
{
  SetName("Fitter") ;
  Z_LOG_NOTE("Fitter constructor called for: " + PtrStr()) ;
}

//--------------------------------------------------------------
// Destructor
Fitter::~Fitter() 
{
  Z_LOG_NOTE("Fitter destructor called for: " + PtrStr()) ;
}

//--------------------------------------------------------------
double Fitter::FitIdx(const size_t& idx) const
{
  //......................................................
  // Data
  TH1F hh("Events", "Events", obs_set[idx].size(), 0, 2370) ;

  // TF1 dataFlat("dataFlat","1",-5,5);
  // hh.FillRandom("dataFlat", 1000) ;

  for(size_t i=0; i<obs_set[idx].size(); i++)
	{
    hh.SetBinContent(i+1, obs_set[idx][i].val) ;
    hh.SetBinError(i+1,  obs_set[idx][i].err) ;
  }

  RooRealVar x("x", "x", 0, 2370) ;
  RooDataHist data("data", "Dataset with x", x, &hh) ;
  //......................................................

  //......................................................
  // Gauss Histogram PDF Model
  // TH1F A_th1f_pdf("A_pdf_th1", "A pdf_th1", 50, -5, 5) ;
  TH1F B_th1f_pdf("B_pdf_th1", "B pdf_th1", obs_set[idx].size(), 0, 2370) ;

  // TF1 A_pdfGauss("A_pdfGauss","TMath::Gaus(x,1,2)",-5,5);
  // TF1 B_pdfFlat("B_pdfFlat","1",-5,5);

  // A_th1f_pdf.FillRandom("A_pdfGauss", 1000) ;
  // B_th1f_pdf.FillRandom("B_pdfFlat", 1000) ;
  for(size_t i=0; i<obs_set[idx].size(); i++)
  {
    B_th1f_pdf.SetBinContent(i+1,1) ;
    // hh.SetBinError(i+1,  bins[i].GetNCount().err) ;
  }

  RooDataHist A_pdf_data("A_pdf_data", "A pdf data", x, &sig_shape_hist[idx]) ;
  RooDataHist B_pdf_data("B_pdf_data", "B pdf data", x, &B_th1f_pdf) ;


  RooHistPdf A_hpdf("A_hist_pdf", "A Histogram pdf", x, A_pdf_data) ;
  RooHistPdf B_hpdf("B_hist_pdf", "B Histogram pdf", x, B_pdf_data) ;
  //......................................................


  //......................................................
  // Composite Model:  f*A(x) + (1-f)*B(x)
  // RooRealVar f("f","A fraction",0.3,0.,1.) ;

  // Extended Likelihood
  RooRealVar A_frac("a_frac","A fraction",10,0.,5000.) ; 
  RooRealVar B_frac("b_frac","B fraction",100,0.,10000.) ;

  RooAddPdf AB_hpdf("Comb_H_pdf", "Combined Hist Pdf", 
                      RooArgList(A_hpdf, B_hpdf), RooArgList(A_frac, B_frac));
  //......................................................

  //......................................................
  // Fitting to data
  // AB_hpdf.fitTo(data) ;
  RooAbsReal *nll = AB_hpdf.createNLL(data);

  // Create MINUIT interface object
  RooMinimizer m(*nll);

  // Activate verbose logging of MINUIT parameter space stepping
  // m.setVerbose(kTRUE);

  // Call MIGRAD to minimize the likelihood
  m.migrad();

  // Print values of all parameters, that reflect values (and error estimates)
  // that are back propagated from MINUIT
  AB_hpdf.getParameters(x)->Print("s");

  // Disable verbose logging
  // m.setVerbose(kFALSE);

  // Run HESSE to calculate errors from d2L/dp2
  m.hesse();

  // Print value (and error) of A_frac parameter, that reflects
  // value and error back propagated from MINUIT
  A_frac.Print();

  // 0.5 for 1-sigma (in 1D)
  m.setErrorLevel(6.1803/2.0) ;

  // Run MINOS on A_frac parameter only
  m.minos(RooArgSet(A_frac, B_frac));

  // Print value (and error) of A_frac parameter, that reflects
  // value and error back propagated from MINUIT
  A_frac.Print();

  //......................................................

  //......................................................
  // S a v i n g   r e s u l t s ,   c o n t o u r   p l o t s
  //......................................................
  // Save a snapshot of the fit result. This object contains the initial
  // fit parameters, the final fit parameters, the complete correlation
  // matrix, the EDM, the minimized FCN , the last MINUIT status code and
  // the number of times the RooFit function object has indicated evaluation
  // problems (e.g. zero probabilities during likelihood evaluation)
  RooFitResult *r = m.save();

#if MOD_MAKE_CONTOUR

  m.setErrorLevel(0.5) ;

  // Run MINOS on A_frac parameter only
  m.minos(RooArgSet(A_frac, B_frac));

  // Make contour plot of mx vs sx at 1,2,3 sigma
  RooPlot *frame = m.contour(A_frac, B_frac, 1, 2, 3);
  frame->SetTitle("Minuit contour plot (Modulation)");

  frame->SetAxisRange(0, 150, "X") ;
  frame->SetAxisRange(800, 1100, "Y") ;
  
  // Plotting the contours
  TCanvas c2("c2", "Contour Canvas", 1000, 800) ;
  frame->Draw() ;
  c2.SaveAs((wrk_dir + "/Modulation_Contours.pdf").Str().c_str());
#endif

  // Print the fit result snapshot
  r->Print("v");
  //......................................................

  //......................................................
  // Fit Results
  double CL_2sigma_A_frac = A_frac.getAsymErrorHi() + A_frac.getVal();

  std::cout << "\n.....................................................\n" ;
  // std::cout << " Upper limit on A_frac = " << CL_2sigma_A_frac << "\n" ;
  std::cout << " 2 sigma C.L. on A_frac["<<idx<<"] = " << CL_2sigma_A_frac  << "\n" ;
  std::cout << " B_frac best fit = " << B_frac.getVal() / obs_set[idx].size() << "\n" ;
  std::cout << ".....................................................\n" ;

  std::cout << " Norm of Signal Shape TH1F = " << sig_shape_hist[idx].Integral() << "\n" ;
  //......................................................

  //......................................................
  // Plotting
  // !! Changing A_frac to it's upper limit for the plots !!
  A_frac.setVal(CL_2sigma_A_frac) ;
  TCanvas c1("c1", "HistPdf Canvas", 1600, 800) ;
  RooPlot* xframe = x.frame() ;
  
  xframe->SetTitle("Constant Observed Events") ;
  xframe->SetXTitle("t (Days)") ;

  data.plotOn(xframe, RooFit::LineColor(kBlue) ) ;
  AB_hpdf.plotOn(xframe, RooFit::LineColor(kPink) ) ;
  AB_hpdf.plotOn(xframe, RooFit::Components("A_hist_pdf"), RooFit::LineColor(kGreen) ) ;
  // AB_hpdf.plotOn(xframe, RooFit::Components("B_hist_pdf"), RooFit::LineStyle(kDotted),
                  // RooFit::LineColor(kCyan) ) ;

  xframe->SetAxisRange(obs_set[idx][0].val - 1.1*obs_set[idx][0].err, 
                        obs_set[idx][0].val + 1.1*obs_set[idx][0].err, "Y") ;

  xframe->Draw() ;

  char out_file_char[150] ;    
  sprintf(out_file_char, "RooFit/Modulation_Hist_%lu.pdf", idx) ;
  
  c1.SaveAs((wrk_dir + "/" + out_file_char).Str().c_str()) ;
  //......................................................

  return CL_2sigma_A_frac / sig_shape_hist[idx].Integral() ;
}

//--------------------------------------------------------------
void Fitter::AddObsCounts(const Bin& b) 
{
  obs_set.push_back(b.GetTBinObsSet()) ;
}

//--------------------------------------------------------------
void Fitter::AddSigShape(const TH1F& sig_shape) 
{
  sig_shape_hist.push_back(sig_shape) ;
}

//--------------------------------------------------------------
std::vector<double> Fitter::Fit()
{
  fit_results.reserve( obs_set.size() ) ;
  for(size_t i=0 ; i < obs_set.size() ; i++)
  {
    fit_results.emplace_back(FitIdx(i)) ;
  }

  return fit_results ;
}
//--------------------------------------------------------------
double Fitter::Fit2D() const
{
  std::vector<Zaki::Math::Quantity> obs_list_1d ;
  obs_list_1d.reserve(obs_set.size()*200) ;

  for (size_t i = 0; i < obs_set.size() ; i++)
  {
    obs_list_1d.insert(obs_list_1d.end(), obs_set[i].begin(), obs_set[i].end());
  }
  
  //......................................................
  // 2D Data
  TH1F hh("Events(AllBins)", "Events(All Bins)", obs_list_1d.size(), 0, 2370*obs_set.size()) ;
  // TF1 dataFlat("dataFlat","1",-5,5);
  // hh.FillRandom("dataFlat", 1000) ;

  // for(size_t i=0; i<obs_set[idx].size(); i++)
	// {
  //   hh.SetBinContent(i+1, obs_set[idx][i].val) ;
  //   hh.SetBinError(i+1,  obs_set[idx][i].err) ;
  // }

  // RooRealVar x("x", "x", 0, 2370) ;
  // RooDataHist data("data", "Dataset with x", x, &hh) ;
  //......................................................
  return 0 ;
}

//==============================================================
