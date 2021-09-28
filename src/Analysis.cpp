/*
  Analysis class

*/

// Creating directory
#include <sys/stat.h>

#include <Math/QuantFuncMathMore.h>
#include <root6/Minuit2/MnScan.h>
#include <root6/Minuit2/MnContours.h>


#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <Zaki/Math/GSLFuncWrapper.hpp>
#include <Zaki/Math/GSLfdfWrapper.hpp>

#include <Zaki/File/CSVIterator.hpp>

// Local headers
// #include "DMSS/ChiSqrd.hpp"
#include "DMSS/LogLikeli.hpp"

#include "DMSS/Analysis.hpp"

//==============================================================

//--------------------------------------------------------------
// Constructor
Analysis::Analysis() 
  : Prog("Analysis", true)
{
  Z_LOG_NOTE("Analysis constructor called for: " + PtrStr()) ;
  m_SatBundles.reserve(10) ;
}

//--------------------------------------------------------------
// Copy Constructor
Analysis::Analysis(const Analysis& other) :
Prog("Analysis", true),
set_satPtr_flag(other.set_satPtr_flag), cpy_ctr_called(other.cpy_ctr_called),
assign_op_called(other.assign_op_called),
sig_shape_hist(other.sig_shape_hist), 
// all_bins(other.all_bins),
focus_bins(other.focus_bins), thresh_limits(other.thresh_limits),
thresh_limits_true(other.thresh_limits_true), e_bins(other.e_bins),
// all_bin_periods(other.all_bin_periods), 
focus_bin_periods(other.focus_bin_periods),
boost_fit_results(other.boost_fit_results), boost_fit_val(other.boost_fit_val),
e_cut_val(other.e_cut_val),
m_SatBundles(other.m_SatBundles)
{
  Z_LOG_NOTE("Analysis copy constructor called: from " + other.PtrStr() + " --> " + PtrStr()) ;
  
  //.............................
  // Copying satellites
  
  if(other.modelPtr != nullptr)
    modelPtr = other.modelPtr->Clone() ;
  //.............................

  //.............................
  // Prog base-class members
  wrk_dir = other.wrk_dir ;
  name        = other.name ;
  set_name_flag = other.set_name_flag    ;
  set_wrk_dir_flag = other.set_wrk_dir_flag ;
  //.............................

  cpy_ctr_called = true ;
}

//--------------------------------------------------------------
// Assignment operator
Analysis& Analysis::operator=(const Analysis& other) 
{

  Z_LOG_NOTE("Analysis '=' operator called: from " + other.PtrStr() + " --> " + PtrStr()) ;

  if(this == &other) return *this ;

  else
  {
    set_satPtr_flag = other.set_satPtr_flag;
    cpy_ctr_called  = other.cpy_ctr_called ;
    assign_op_called  = other.assign_op_called ;
    sig_shape_hist= other.sig_shape_hist;
    focus_bins= other.focus_bins;
    thresh_limits= other.thresh_limits;
    thresh_limits_true= other.thresh_limits_true; e_bins= other.e_bins;
    focus_bin_periods= other.focus_bin_periods;
    boost_fit_results= other.boost_fit_results; 
    boost_fit_val= other.boost_fit_val;
    e_cut_val= other.e_cut_val;
    m_SatBundles = other.m_SatBundles ;
    // Pointer member variables
    modelPtr = other.modelPtr->Clone() ;
    assign_op_called = true ;

    return *this ;
  }

}

//--------------------------------------------------------------
// Destructor
Analysis::~Analysis() 
{
  Z_LOG_NOTE("Analysis destructor called for: " + PtrStr()) ;
}

//--------------------------------------------------------------
// Overriding the base method
void Analysis::SetWrkDir(const Zaki::String::Directory& input) 
{
  wrk_dir = input ;
  set_wrk_dir_flag = true ;

  // ............ Creating a directory ............
  if (mkdir(input.Str().c_str(), ACCESSPERMS) == -1) 
  {
    Z_LOG_INFO("Directory '"+input.Str()+"' wasn't created, because: "+strerror(errno)+".") ;  
  }
  else
    Z_LOG_INFO("Directory '" + wrk_dir.Str() + "' created."); 
  // .................................................

  for(auto satB : m_SatBundles)
  {
    satB->SetWrkDir(input) ;
  }
}

//--------------------------------------------------------------
void Analysis::AttachSatellite(std::shared_ptr<Satellite> s)
{
  
  if(set_wrk_dir_flag)
    s->SetWrkDir(wrk_dir.Str().c_str()) ;

  m_SatBundles.emplace_back(s) ;
  
  set_satPtr_flag = true ;
}

//--------------------------------------------------------------
std::shared_ptr<Satellite> Analysis::GetSatellite(const size_t& i)
{
  if(m_SatBundles.size() - 1 < i )
  {
    Z_LOG_ERROR("Satellite index '"+ std::to_string(i) + "' is out of bound!") ;
    return nullptr;
  }
  return m_SatBundles[i].GetSatellite() ;
}

//--------------------------------------------------------------
void Analysis::SetECut(double energy)
{
  e_cut_val = energy;

  // If model is already added
  if (modelPtr)
    modelPtr->SetECut(energy); 
}

//--------------------------------------------------------------
double Analysis::GetECut() const
{
  return e_cut_val;
}

//--------------------------------------------------------------
// Threshold fitting (without minuit)
void Analysis::DoThresh(const size_t& sat_idx, 
  std::vector<size_t>* in_bins_idx, Zaki::Math::Grid2D* grid_in)
{

  using namespace CONFIND ;

  if(GetSatellite(sat_idx)->GetData().size() == 0)
  {
    Z_LOG_ERROR("Satellite data is empty!") ;
    exit(EXIT_FAILURE) ;
  }
  
  e_bins.clear() ;
  //.........................................
  std::vector<std::vector<double>> t_mu95 ; t_mu95.reserve(GetBinPeriod(sat_idx).size()+1) ;
  std::vector<double> mu95 ;  mu95.reserve(GetBinPeriod(sat_idx).size()+1) ;
  double tmp_mu95 ;

  if(in_bins_idx)
  {
    size_t j = 0 ;
    for(size_t i : *in_bins_idx)
    {
      t_mu95.push_back({(double)GetBinPeriod(sat_idx)[j], tmp_mu95}) ;
      mu95.push_back(tmp_mu95) ;
      e_bins.push_back(m_SatBundles[sat_idx].GetTimeBin(i,j).GetECenter().val);
      j++ ;
    }
  } 
  else
  {
    for(size_t i=0 ; i < m_SatBundles[sat_idx]->GetData().GetBins().size() ; ++i)
    {
      t_mu95.push_back({(double)GetBinPeriod(sat_idx)[i], tmp_mu95}) ;
      mu95.push_back(tmp_mu95) ;
      e_bins.push_back(m_SatBundles[sat_idx].GetTimeBin(i,i).GetECenter().val);
    }
  }


  Zaki::File::VecSaver my_saver(m_SatBundles[sat_idx]->GetWrkDir() + "/LikeLi/t_mu95_" + m_SatBundles[sat_idx]->GetName() + ".dat", 
                                Zaki::File::FileMode::Write) ;
  my_saver.Export2D(t_mu95) ;
  //.........................................

  // Initializing the model
  if(modelPtr)
    modelPtr->Init() ;
  else
  {
    Z_LOG_ERROR("Model is not set, use 'SetModel' first.") ;
    return;
  }
  
  // ...........................
  // Setting up the grid
  // ...........................
  MemFuncContWrapper<Analysis, double (Analysis::*) (double, double)>* mfcwPtr;

  switch (plot_mode)
  {
    case PlotMode::Ldec_Gann:
      mfcwPtr = new MemFuncContWrapper<Analysis,
                double (Analysis::*) (double, double)>(*this, &Analysis::ContFuncThresh_LG) ;
      break;
    
    case PlotMode::Mdm_Gann:
      mfcwPtr = new MemFuncContWrapper<Analysis,
                double (Analysis::*) (double, double)>(*this, &Analysis::ContFuncThresh_MG) ;
      break;

    case PlotMode::Ldec_Mdm:
      mfcwPtr = new MemFuncContWrapper<Analysis, 
                double (Analysis::*) (double, double)>(*this, &Analysis::ContFuncThresh_LM) ;
      break;

    default:
      break;
  }

  //..................
  // Checking the grid input
  if (grid_in)
  {
    Z_LOG_INFO("Using the grid input from 'Analysis'.") ;
    (*mfcwPtr)->SetGrid(*grid_in) ;
  }
  // Checking the grid in the model
  else if(modelPtr->GetTGrid(plot_mode))
  {
    Z_LOG_INFO("Using the grid from 'Model'.") ;
    (*mfcwPtr)->SetGrid(modelPtr->GetTGrid(plot_mode)->grid) ;
  }
  // Emitting errors and stopping the process
  else
  {
    Z_LOG_ERROR("Analysis aborted: No grid was input, and there is also none in the model!") ;
    return ;
  }
  //..................
  std::string tmp_name = "Combined" ; 
  if(m_SatBundles.size() == 1)
    tmp_name = m_SatBundles[sat_idx]->GetName() ;

  (*mfcwPtr)->SetWrkDir(wrk_dir + "/" + tmp_name + "/Fit") ;

  (*mfcwPtr)->SetThreads(8) ;

  // ...........................
  // Saving time by checking if the grid values are the
  //  same for various contours
  if(modelPtr->GetTGrid(plot_mode)->fixed)
  {
    std::vector<std::string> labels ;
    for (size_t i = 0; i < e_bins.size(); ++i)
    { 
      if(in_bins_idx)
        labels.emplace_back("T_"+std::to_string(GetBinPeriod(sat_idx)[i])) ;
      else
        labels.emplace_back("M_" +std::to_string((int)e_bins[i]*modelPtr->GetDecayProd().size())) ; 
    }

    (*mfcwPtr)->SetContVal(mu95, labels) ;
    (*mfcwPtr)->SetGridVals(CONFIND::ContourFinder::Mode::Fast) ;
  }
  else
  {
    for (size_t i = 0; i < mu95.size(); ++i)
    {
      switch (plot_mode)
      {
        case PlotMode::Ldec_Gann:
        // Changing the DM mass in the model
        modelPtr->SetDMMass(e_bins[i]*modelPtr->GetDecayProd().size()) ;
        break;

        case PlotMode::Mdm_Gann:

        break;

        case PlotMode::Ldec_Mdm:

        break;

        default:
        break;
      }

      if(in_bins_idx)
      {
        modelPtr->SetActiveBin(m_SatBundles[sat_idx]->GetData().GetBins()[(*in_bins_idx)[i]]) ;
        (*mfcwPtr)->SetContVal({mu95[i]}, {"T_"+std::to_string(GetBinPeriod(sat_idx)[i])}) ;
      }
      else
      {
        modelPtr->SetActiveBin(m_SatBundles[sat_idx]->GetData().GetBins()[i]) ;
        (*mfcwPtr)->SetContVal({mu95[i]}, {"M_"+std::to_string((int)e_bins[i]*modelPtr->GetDecayProd().size())}) ;
      }

      // Need to update mfw after any change
      UpdateMFCW(mfcwPtr) ;

      (*mfcwPtr)->SetGridVals(CONFIND::ContourFinder::Mode::Fast) ;
    }
  }
  // ...........................

  // ...........................
  // Generating results
  // ...........................
  // ............ Creating a directory ............
  if (mkdir((wrk_dir + "/" + tmp_name + "/Fit/Contour Points/Threshold").Str().c_str(), ACCESSPERMS) == -1) 
  {
    Z_LOG_INFO("Directory '"+ (wrk_dir + "/" + tmp_name + "/Fit/Contour Points/Threshold").Str()+"' wasn't created, because: "+strerror(errno)+".") ;    
  }
  else
    Z_LOG_INFO("Directory '"+(wrk_dir + "/" + tmp_name + "/Fit/Contour Points/Threshold").Str()+"' created."); 
  // .................................................
  (*mfcwPtr)->ExportContour("Contour Points/Threshold/"+ modelPtr->GetName() + "_threshold_cont",  Zaki::File::FileMode::Write) ; 

  size_t tmp_bins_size = 0 ; 
  double tmp_leg_x     = 0 ;
  std::string tmp_leg_label = "M_{DM} GeV";
  if(in_bins_idx) 
  { 
    tmp_bins_size = in_bins_idx->size() ;
    tmp_leg_x = 0.05 ;
    tmp_leg_label = "T (Day)" ;
  } 
  else            
  { 
    tmp_bins_size = m_SatBundles[sat_idx]->GetData().GetBins().size() ; 
  }
  // ...........................
  // Making the plot legend
  (*mfcwPtr)->MakeLegend(true, tmp_leg_label.c_str(), "user") ;
  (*mfcwPtr)->GetLegend()->SetX1(0.75 + tmp_leg_x) ; (*mfcwPtr)->GetLegend()->SetY1(0.10) ;
  (*mfcwPtr)->GetLegend()->SetX2(0.90) ; (*mfcwPtr)->GetLegend()->SetY2(0.15 + tmp_bins_size*0.03) ;
  (*mfcwPtr)->GetLegend()->SetTextSize(0.025) ;

  std::vector<std::string> tmp_labels ;
  char tmp[200] ;

  for (size_t i = 0; i < tmp_bins_size; i++)
  {
    if(in_bins_idx)
      sprintf(tmp, " %d", m_SatBundles[sat_idx].GetBinPeriod()[i]) ;
    else
      sprintf(tmp, "%.1e", e_bins[i]*modelPtr->GetDecayProd().size()) ;

    tmp_labels.emplace_back(tmp) ;
  }
    
  
  (*mfcwPtr)->SetLegendLabels(tmp_labels) ;
  // ...........................

  if(modelPtr->GetTGrid(plot_mode)->option.joined)
    (*mfcwPtr)->SetPlotConnected() ;
  
  (*mfcwPtr)->SetWidth(modelPtr->GetTGrid(plot_mode)->option.w) ;
  (*mfcwPtr)->SetHeight(modelPtr->GetTGrid(plot_mode)->option.h) ;

  sprintf(tmp,"%s - Threshold (%s)", modelPtr->GetName().c_str(),
          m_SatBundles[sat_idx]->GetName().c_str()) ;
  (*mfcwPtr)->Plot(modelPtr->GetName() + "_Threshold", tmp, 
           modelPtr->GetTGrid(plot_mode)->option.xLabel.c_str(), 
           modelPtr->GetTGrid(plot_mode)->option.yLabel.c_str()) ;
  // ...........................

  delete mfcwPtr;
  // ................................................
}

//--------------------------------------------------------------
/// Scans the likelihood function vs parameters
std::vector<double> Analysis::ScanParThresh(const size_t& b_idx, const size_t& t_idx, const size_t& sat_idx) 
{
  Bin b = m_SatBundles[sat_idx].GetTimeBin(b_idx,t_idx) ;

  TH1F sig_shape = m_SatBundles[sat_idx]->GetSigShape(b.GetECenter().val) ;

  ROOT::Minuit2::LogLikeli  fcn;

  fcn.AddObsCounts(b) ;
  fcn.AddSigShape(sig_shape)  ;

  using namespace CONFIND ;
  MemFuncContWrapper<ROOT::Minuit2::LogLikeli, 
                      double (ROOT::Minuit2::LogLikeli::*) (double, double)> 
                      mfcw(fcn, &ROOT::Minuit2::LogLikeli::ThreshCont) ;


    //  ( PlotMode::Ldec_Gann,
    //        {{{5e4, 1e12}, tmp_Res, "Log"}, {{3e-6, 9e-1}, tmp_Res, "Log"}},
    //        {1000, 1000, "L_{dec} [km]", "#Gamma [GeV]", true},
    //       true)
  double tmp_bg = b.GetTBinObsSet()[0].val  ;

//  mfcw->SetGrid({{{1e-5, 3e-2}, 400, "Log"}, {{0 , tmp_bg * 1.1 }, 400, "Linear"}}) ;
 mfcw->SetGrid({{{1e-5, 3e-2}, 400, "Log"}, {{tmp_bg*0.1 , tmp_bg * 1.1 }, 400, "Linear"}}) ;

  
  mfcw->SetWrkDir(wrk_dir + "/" + m_SatBundles[sat_idx]->GetName()  + "/LikeLi") ;
  
  mfcw->SetContVal({ROOT::MathMore::chisquared_quantile(0.68, 2),
                    ROOT::MathMore::chisquared_quantile(0.90, 2),
                    ROOT::MathMore::chisquared_quantile(0.95, 2)}) ;

  mfcw->SetGridVals(CONFIND::ContourFinder::Mode::Fast) ;

  mfcw->MakeLegend(true, "#Delta [ -2Log(L) ]", "user") ;
  mfcw->GetLegend()->SetX1(0.10) ; mfcw->GetLegend()->SetY1(0.10) ;
  mfcw->GetLegend()->SetX2(0.25) ; mfcw->GetLegend()->SetY2(0.30) ;
  mfcw->GetLegend()->SetTextSize(0.025) ;

  mfcw->Plot("Thresh/Plots/Threshold_B" + std::to_string(b_idx) + "_T="
             + std::to_string(GetBinPeriod(sat_idx)[t_idx]), ("-2 Log Likelihood ( T = "
             + std::to_string(GetBinPeriod(sat_idx)[t_idx]) +" )").c_str(),
            "#mu",
            "bg") ;
  mfcw->ExportContour("Thresh/Conts/B=" + std::to_string(b_idx)  
                      + "_T=" + std::to_string(GetBinPeriod(sat_idx)[t_idx]), 
                      Zaki::File::FileMode::Write ) ;

  // Removing the contours we don't need
  for (auto &&i : {"2.28e+00", "4.61e+00"})
  {
    std::remove((wrk_dir + "/" + m_SatBundles[sat_idx]->GetName() 
                                    + "/LikeLi/Thresh/Conts/B=" + std::to_string(b_idx)  
                  + "_T=" + std::to_string(GetBinPeriod(sat_idx)[t_idx]) + "_" + i).Str().c_str()) ;
  }



  std::ifstream file( (wrk_dir + "/" + m_SatBundles[sat_idx]->GetName() 
                      + "/LikeLi/Thresh/Conts/B=" + std::to_string(b_idx)  
                      + "_T=" + std::to_string(GetBinPeriod(sat_idx)[t_idx]) + "_5.99e+00").Str());

  // Error opening the file
  if (file.fail()) 
  {
    Z_LOG_ERROR("File '"+(wrk_dir + "/" + m_SatBundles[sat_idx]->GetName() 
                                    + "/LikeLi/"
                                    + "Thresh/Conts/B=" + std::to_string(b_idx)  
                      + "_T=" + std::to_string(GetBinPeriod(sat_idx)[t_idx]) + "_5.99e+00").Str() +"' cannot be opened!") ;
    exit(EXIT_FAILURE) ;
    return {-1, -1, -1};
  }

  std::vector<double> mu95_vals ; mu95_vals.reserve(500) ;

  // Reading the input file
  for(Zaki::File::CSVIterator loop(file); loop != Zaki::File::CSVIterator(); ++loop)
  {
    mu95_vals.emplace_back(std::stof((*loop)[0])) ;
  }

return {(double)b_idx, (double)GetBinPeriod(sat_idx)[t_idx], *std::max_element(mu95_vals.begin(), mu95_vals.end())} ;

}
//--------------------------------------------------------------
// Plots the parameter space given the mu95 input file
void Analysis::PlotParamSpaceThresh(const size_t& sat_idx, const bool bypass) 
{
  std::vector<std::vector<double>> tmp_vec ;
 if(!bypass)
 {
    for (size_t i = 0; i < GetBinPeriod(0).size(); i++)
    {
      tmp_vec.emplace_back( ScanParThresh(i, i, 0) ) ;
    }

    Zaki::File::VecSaver my_saver(wrk_dir + "/" + m_SatBundles[sat_idx]->GetName() 
                                      + "/LikeLi/"
                                      + "Thresh/Conts/"+ m_SatBundles[sat_idx]->GetName() 
                                      +"_Thresh.dat", 
                                  Zaki::File::FileMode::Write) ;
    my_saver.Export2D(tmp_vec) ;
  }
  // Limits are already found we just need to read them from file
  else
  {
    std::ifstream file( (wrk_dir + "/" + m_SatBundles[sat_idx]->GetName() 
                                      + "/LikeLi/"
                                      + "Thresh/Conts/"+ m_SatBundles[sat_idx]->GetName() 
                                      +"_Thresh.dat" ).Str());

    // Error opening the file
    if (file.fail()) 
    {
      Z_LOG_ERROR("File '"+(wrk_dir + "/" + m_SatBundles[sat_idx]->GetName() 
                                      + "/LikeLi/"
                                      + "Thresh/Conts/"+ m_SatBundles[sat_idx]->GetName() 
                                      +"_Thresh.dat" ).Str() +"' cannot be opened!") ;
      exit(EXIT_FAILURE) ;
      return ;
    }

    // Reading the input file
    for(Zaki::File::CSVIterator loop(file); loop != Zaki::File::CSVIterator(); ++loop)
    {
      tmp_vec.push_back({std::stof((*loop)[0]), std::stof((*loop)[1]), std::stof((*loop)[2])}) ;
    }
  }
  
  //....................................
  // Initializing the model
  if(modelPtr)
    modelPtr->Init() ;
  else
  {
    Z_LOG_ERROR("Model is not set, use 'SetModel' first.") ;
    return;
  }
  //....................................

  using namespace CONFIND ;
  MemFuncContWrapper<Analysis, 
                      double (Analysis::*) (double, double)> 
                      mfcw(*this, &Analysis::ContFuncThresh_LG) ;


  mfcw->SetGrid(modelPtr->GetTGrid(plot_mode)->grid) ;
  
  mfcw->SetWrkDir(wrk_dir + "/" + m_SatBundles[sat_idx]->GetName()  + "/LikeLi") ;
    

  std::vector<double> cont_vals ; cont_vals.reserve(20) ;
  std::vector<std::string> cont_labels ; cont_labels.reserve(20) ;

  // cont_labels = {"200", "230", "270", "310", "360", "410", "470", "540", "620",
  // "710", "810",  "930", "1100", "1200", "1400", "1600", "1900", "2100",
  // "2500", "2800", "3200" , "3700", "4300", "4900", "5600", "6500",
  // "7400", "8500"} ;
  char tmp_lab[100] ;
  for (size_t i = 0; i < tmp_vec.size(); i++)
  {
    cont_vals.push_back(tmp_vec[i][2]) ;
    sprintf(tmp_lab, "%d", 10*(int)((m_SatBundles[sat_idx].GetTimeBin(i,i).GetECenter().val*modelPtr->GetDecayProd().size())/10)) ;

    // Changing the DM mass in the model
    modelPtr->SetDMMass(m_SatBundles[sat_idx].GetTimeBin(i,i).GetECenter().val*modelPtr->GetDecayProd().size()) ;
    modelPtr->SetActiveBin(m_SatBundles[sat_idx]->GetData().GetBins()[i]) ;
    mfcw->SetContVal({tmp_vec[i][2]}, {tmp_lab}) ;
    UpdateMFCW(&mfcw) ;
    mfcw->SetGridVals(CONFIND::ContourFinder::Mode::Fast) ;

    cont_labels.push_back(tmp_lab) ;
  }




  mfcw->MakeLegend(true, "M_{DM}  [GeV]", "user") ;
  mfcw->GetLegend()->SetX1(0.75) ; mfcw->GetLegend()->SetY1(0.10) ;
  mfcw->GetLegend()->SetX2(0.90) ; mfcw->GetLegend()->SetY2(0.15 + tmp_vec.size()*0.03) ;
  mfcw->GetLegend()->SetTextSize(0.025) ;
  mfcw->SetLegendLabels(cont_labels) ;
  
  char tmp_plot_label[200] ;
  sprintf(tmp_plot_label,"%s - Threshold (%s)", 
          modelPtr->GetName().c_str(), m_SatBundles[sat_idx]->GetName().c_str()) ;  

  mfcw->Plot(m_SatBundles[sat_idx]->GetName() + "_Threshold_" 
            //+ "_B=0_TCompar", 
             + modelPtr->GetName(),
              tmp_plot_label,
              modelPtr->GetTGrid(plot_mode)->option.xLabel.c_str(), 
              modelPtr->GetTGrid(plot_mode)->option.yLabel.c_str()) ;

  mfcw->ExportContour("Thresh/Conts/" + modelPtr->GetName()+ "/" + 
                      m_SatBundles[sat_idx]->GetName() + "_" + modelPtr->GetName() 
                      // +"_B=0" // + std::to_string(b_idx)  
                      // + "_TCompar"
                      , 
                      Zaki::File::FileMode::Write ) ;
}

//--------------------------------------------------------------
/// Scans the likelihood function vs parameters (done by Minuit)
void Analysis::ScanParMinuit(const size_t& b_idx, const size_t& sat_idx) 
{
  Bin b = m_SatBundles[sat_idx].GetTimeBin(b_idx,b_idx) ;

  TH1F sig_shape = m_SatBundles[sat_idx]->GetSigShape(b.GetECenter().val) ;

  ROOT::Minuit2::LogLikeli  fcn;

  fcn.AddObsCounts(b) ;
  fcn.AddSigShape(sig_shape)  ;

  // ****************
  //  Minimize
  // ****************
  //(para_name, para_start_value, para_uncertainties)
  ROOT::Minuit2::MnUserParameters upar ;
  upar.Add("mu", 0.01, 1e-10);
  upar.Add("bg", b.GetTBinObsSet()[0].val, 1, 0, 500);
  // upar.Add("period", bin_per, 0.00000001);

  // set limits
  upar.SetLowerLimit("mu", 0);


  // Creat MIGRAD minimizer
  ROOT::Minuit2::MnMigrad migrad(fcn, upar, 2);

  // migrad.Fix("mu") ;

  // // Minimize
  // ROOT::Minuit2::FunctionMinimum min = migrad();

  // std::cout<<"minimum: "<<min<<std::endl;

  // migrad.Release("mu") ;

  std::cout<<"bg: "<<upar.Parameter(1).Value() <<std::endl;


  ROOT::Minuit2::MnScan scan(fcn, upar) ;
  std::vector<std::pair<double, double>> points = scan.Scan(0, 100, 1e-10, 1e-8) ;

  for (auto &&i : points)
  {
    printf("- %5.3e\t %.18e\n", i.first, i.second/(-7.78377063e+04) - 1) ;
  }
  
}

//--------------------------------------------------------------
/// Scans the likelihood function vs parameters
std::vector<double> Analysis::ScanParBoost(const double& in_mass, const size_t& sat_idx) 
{

  Zaki::Math::Range<double> e_range = {GetECut(), in_mass} ; 

  std::vector<double> tmp_best_fit = {1e-7} ;
  ROOT::Minuit2::LogLikeli  fcn ;
  std::vector<double> tmp_out ;

  std::vector<double> fixed_bg_set  ;
  

  for(size_t i=0 ; i < m_SatBundles[sat_idx]->GetData().GetBins().size() ; ++i)
  {
    Bin b = m_SatBundles[sat_idx].GetTimeBin(i,i) ;
    // ..............................
    // out of spectrum condition
    if( b.GetERange().min >= e_range.max || b.GetERange().max <= e_range.min )
      continue ; 


    TH1F sig_shape = m_SatBundles[sat_idx]->GetSigShape(b.GetECenter().val) ;

    // Scaling the histogram by the energy-dependent 
    //  spectrum of the model
    sig_shape.Scale(modelPtr->GetSpectrum(b)) ;
          
    fcn.AddObsCounts(b) ;
    fcn.AddSigShape(sig_shape)  ;
    tmp_best_fit.push_back(b.GetTBinObsSet()[0].val) ; // obs 
    fixed_bg_set.push_back(b.GetTBinObsSet()[0].val) ;
  }

  for (size_t ign_idx = 0; ign_idx < fixed_bg_set.size() ; ign_idx++)
  {
    fcn.SetNDBgSet(ign_idx, fixed_bg_set) ;
    fcn.SetNDBestFit(tmp_best_fit) ;

    using namespace CONFIND ;
    MemFuncContWrapper<ROOT::Minuit2::LogLikeli, 
                        double (ROOT::Minuit2::LogLikeli::*) (double, double)> 
                        mfcw(fcn, &ROOT::Minuit2::LogLikeli::BoostCont) ;

    mfcw->SetGrid({{{1e-5, 5e-3}, 200, "Linear"}, {{fixed_bg_set[ign_idx]*0 , 11 
    //fixed_bg_set[ign_idx] * 1.2 
    }, 200, "Linear"}}) ;
    // mfcw->SetGrid({{{0, 3.5e-3}, 200, "Linear"}, {{30, 34}, 200, "Linear"}}) ;


    mfcw->SetWrkDir(wrk_dir + "/" + m_SatBundles[sat_idx]->GetName()  + "/LikeLi") ;
    
    mfcw->SetContVal({ROOT::MathMore::chisquared_quantile(0.68, tmp_best_fit.size()),
                      ROOT::MathMore::chisquared_quantile(0.90, tmp_best_fit.size()),
                      ROOT::MathMore::chisquared_quantile(0.95, tmp_best_fit.size())}) ;
    mfcw->SetGridVals(CONFIND::ContourFinder::Mode::Fast) ;

    mfcw->MakeLegend(true, "#Delta [ -2Log(L) ]", "user") ;
    mfcw->GetLegend()->SetX1(0.75) ; mfcw->GetLegend()->SetY1(0.10) ;
    mfcw->GetLegend()->SetX2(0.90) ; mfcw->GetLegend()->SetY2(0.30) ;
    mfcw->GetLegend()->SetTextSize(0.025) ;
    mfcw->GetGraph()->GetXaxis()->SetMaxDigits(2) ;

    char tmp[200] ;
    sprintf(tmp, "Boost/Plots/M_%.0f/Boosted_%zu_M=%.0f", in_mass, ign_idx, in_mass) ;
    mfcw->Plot(tmp,
                ("-2 Log Likelihood ( B = " + std::to_string(ign_idx) + ")").c_str(),
                "#mu", "bg") ;
    sprintf(tmp, "Boost/Conts/M_%.0f/Boosted_%zu_M=%.0f",  in_mass, ign_idx, in_mass) ;
    
    // mfcw->ExportContour(tmp, Zaki::File::FileMode::Write ) ;

    
    // Removing the contours we don't need
    for (const double i : {0.68, 0.90})
    {
      sprintf(tmp, "/LikeLi/Boost/Conts/M_%.0f/Boosted_%zu_M=%.0f_%.2e", 
        in_mass, ign_idx, in_mass,
        ROOT::MathMore::chisquared_quantile(i, tmp_best_fit.size()) ) ;

      std::remove((wrk_dir + "/" + m_SatBundles[sat_idx]->GetName() + std::string(tmp)).Str().c_str()) ;
    }

    sprintf(tmp, "/LikeLi/Boost/Conts/M_%.0f/Boosted_%zu_M=%.0f_%.2e", 
              in_mass, ign_idx, in_mass,
              ROOT::MathMore::chisquared_quantile(0.95, tmp_best_fit.size()) ) ;

    std::ifstream file( (wrk_dir + "/" + m_SatBundles[sat_idx]->GetName() 
                        + std::string(tmp) ).Str());

    // Error opening the file
    if (file.fail()) 
    {
      Z_LOG_ERROR("File '"+(wrk_dir + "/" + m_SatBundles[sat_idx]->GetName() 
                            +  std::string(tmp) ).Str() +"' cannot be opened!") ;
      exit(EXIT_FAILURE) ;
      return {-1, -1};
    }

    std::vector<double> mu95_vals ; mu95_vals.reserve(500) ;

    // Reading the input file
    for(Zaki::File::CSVIterator loop(file); loop != Zaki::File::CSVIterator(); ++loop)
    {
      mu95_vals.emplace_back(std::stof((*loop)[0])) ;
    }

    std::cout << "\n #Bins = " << tmp_best_fit.size() << ",\tmu_95 = "
              << *std::max_element(mu95_vals.begin(), mu95_vals.end()) << "\n" ;
    tmp_out.push_back(*std::max_element(mu95_vals.begin(), mu95_vals.end())) ;
  }

  return { in_mass,  *std::max_element(tmp_out.begin(), tmp_out.end())  } ;

}

//--------------------------------------------------------------
void Analysis::PlotParamSpaceBoost(const size_t& sat_idx, const bool bypass) 
{
  
  //....................................
  // Initializing the model
  if(modelPtr)
    modelPtr->Init() ;
  else
  {
    Z_LOG_ERROR("Model is not set, use 'SetModel' first.") ;
    return;
  }
  //....................................

  std::vector<std::vector<double>> tmp_vec ;
  if(!bypass)
  {
    for (size_t i = 0; i < modelPtr->GetDMMassSet().size(); i++)
    {

      // Changing the DM mass
      modelPtr->SetDMMass(modelPtr->GetDMMassSet()[i]) ;

      tmp_vec.emplace_back( ScanParBoost(modelPtr->GetDMMass(), sat_idx) ) ;
    
    }  
    
    Zaki::File::VecSaver my_saver(wrk_dir + "/" + m_SatBundles[sat_idx]->GetName() 
                                    + "/LikeLi/"
                                    + "Boost/Conts/"+ m_SatBundles[sat_idx]->GetName() 
                                    +"_Boost.dat", 
                                Zaki::File::FileMode::Write) ;
    my_saver.Export2D(tmp_vec) ;
  }
  // Limits are already found we just need to read them from file
  else
  {
    std::ifstream file( (wrk_dir + "/" + m_SatBundles[sat_idx]->GetName() 
                                      + "/LikeLi/"
                                      + "Boost/Conts/"+ m_SatBundles[sat_idx]->GetName() 
                                      +"_Boost.dat" ).Str());

    // Error opening the file
    if (file.fail()) 
    {
      Z_LOG_ERROR("File '"+(wrk_dir + "/" + m_SatBundles[sat_idx]->GetName() 
                                      + "/LikeLi/"
                                      + "Boost/Conts/"+ m_SatBundles[sat_idx]->GetName() 
                                      +"_Boost.dat" ).Str() +"' cannot be opened!") ;
      exit(EXIT_FAILURE) ;
      return ;
    }

    // Reading the input file
    for(Zaki::File::CSVIterator loop(file); loop != Zaki::File::CSVIterator(); ++loop)
    {
      tmp_vec.emplace_back(std::stof((*loop)[0]), std::stof((*loop)[1])) ;
    }
  }

return  ;
  using namespace CONFIND ;
  MemFuncContWrapper<Analysis, 
                      double (Analysis::*) (double, double)> 
                      mfcw(*this, &Analysis::ContFuncBoost_LG) ;


  mfcw->SetGrid(modelPtr->GetBGrid(plot_mode)->grid) ;
  
  mfcw->SetWrkDir(wrk_dir + "/" + m_SatBundles[sat_idx]->GetName()  + "/LikeLi") ;

  std::vector<double> cont_vals ; cont_vals.reserve(20) ;
  std::vector<std::string> cont_labels ; cont_labels.reserve(20) ;

  char tmp_lab[100] ;
  for (size_t i = 0; i < tmp_vec.size(); i++)
  {
    cont_vals.push_back(tmp_vec[i][1]) ;
    sprintf(tmp_lab, "%.1e", modelPtr->GetDMMassSet()[i]) ;

    //   Changing the DM mass in the model
    modelPtr->SetDMMass(modelPtr->GetDMMassSet()[i]) ;

    mfcw->SetContVal({tmp_vec[i][1]}, {tmp_lab}) ;
    UpdateMFCW(&mfcw) ;
    mfcw->SetGridVals(CONFIND::ContourFinder::Mode::Fast) ;

    cont_labels.push_back(tmp_lab) ;
  }

  mfcw->MakeLegend(true, "M_{DM}  [GeV]", "user") ;
  mfcw->GetLegend()->SetX1(0.75) ; mfcw->GetLegend()->SetY1(0.10) ;
  mfcw->GetLegend()->SetX2(0.90) ; mfcw->GetLegend()->SetY2(0.15 + tmp_vec.size()*0.03) ;
  mfcw->GetLegend()->SetTextSize(0.025) ;
  mfcw->SetLegendLabels(cont_labels) ;
  mfcw->GetGraph()->GetXaxis()->SetMaxDigits(2) ;
  
  char tmp_plot_label[200] ;
  sprintf(tmp_plot_label,"%s - Boost (%s)", 
          modelPtr->GetName().c_str(), m_SatBundles[sat_idx]->GetName().c_str()) ;  

  mfcw->Plot(m_SatBundles[sat_idx]->GetName() + "_Boost_" 
             + modelPtr->GetName(),
              tmp_plot_label,
              modelPtr->GetBGrid(plot_mode)->option.xLabel.c_str(), 
              modelPtr->GetBGrid(plot_mode)->option.yLabel.c_str()) ;

  mfcw->ExportContour("Boost/Conts/" + modelPtr->GetName()+ "/" + 
                      m_SatBundles[sat_idx]->GetName() + "_" + modelPtr->GetName(), 
                      Zaki::File::FileMode::Write ) ;
}

//--------------------------------------------------------------
// Fitting the histograms and finding the errors (C.L.)
void Analysis::FitBin(Bin b, const size_t& sat_idx)
{
 
  Z_LOG_INFO("Performing the threshold fit for E = " 
              + std::to_string(b.GetECenter().val) + " GeV.") ;

  TH1F sig_shape = m_SatBundles[sat_idx]->GetSigShape(b.GetECenter().val) ;

  ROOT::Minuit2::LogLikeli   fcn;

  fcn.AddObsCounts(b) ;
  fcn.AddSigShape(sig_shape)  ;

  // ****************
  //  Minimize
  // ****************
  //(para_name, para_start_value, para_uncertainties)
  ROOT::Minuit2::MnUserParameters upar ;
  upar.Add("mu", 3e-4, 1e-7);
  upar.Add("bg", b.GetTBinObsSet()[0].val, 1);

  std::cout << "\n --------> b.GetTBinObsSet()[0].val: " << b.GetTBinObsSet()[0].val << "\n" ;
  // upar.Add("period", bin_per, 0.00000001);

  // set limits
  upar.SetLowerLimit("mu", 0);


  // Creat MIGRAD minimizer
  int strategy_counter = 2 ;
  ROOT::Minuit2::MnMigrad migrad(fcn, upar, strategy_counter);

  // Minimize
  ROOT::Minuit2::FunctionMinimum min = migrad();

  //.................
  // If minimization failed:
  while(!min.IsValid() && strategy_counter < 2) 
  {
    strategy_counter++ ;
    // try with higher strategy
    Z_LOG_ERROR("Migrad min is invalid, trying with strategy = "+ std::to_string(strategy_counter)+ ".");
    ROOT::Minuit2::MnMigrad migrad(fcn, upar, strategy_counter);
    min = migrad();
  }
  //.................

  // output
  std::cout<<"minimum: "<<min<<std::endl;

  // ****************
  //  Error analysis
  // ****************
  // MINOS Error analysis
  ROOT::Minuit2::MnMinos Minos(fcn, min);

  // MINOS errors
  // -2Log & chi2 have the same normalization:
  // Zaki::Math::NDimContLevel cont_level((int)upar.Params().size(), conf_level) ;
  // fcn.SetErrorDef(cont_level.Solve().Up);
  fcn.SetErrorDef(ROOT::MathMore::chisquared_quantile(conf_level,
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

  sprintf(tmp_char, " ---> Minos error for E=%.1f, t = %d is invalid!", b.GetECenter().val, sig_shape.GetSize()) ;
  if(!e0.IsValid()) Z_LOG_ERROR(tmp_char) ;

  double out = Histgen(sat_idx, e0, e1, b, sig_shape) ;

  // thresh_limits_true.push_back(
  //   (e0.Min() + e0.Upper())
  //   *m_SatBundles[sat_idx]->GetNumScaling() / m_SatBundles[sat_idx]->GetExpNorm() 
  //   ) ;

  thresh_limits_true.push_back(
    (e0.Min() + e0.Upper()) * m_SatBundles[sat_idx]->GetNumScaling()  ) ;

  thresh_limits.push_back(out);
  e_bins.push_back(b.GetECenter().val);
    // fitted_bins.push_back(b) ;
}

//--------------------------------------------------------------
std::vector<double> Analysis::GetThreshLimits() const
{
  return thresh_limits_true ;
}

//--------------------------------------------------------------
/// Threshold Case ( bin-by-bin )
/// If in_bins_idx is nullptr, it will fit all the satelite's bins
void Analysis::FitThreshold(const size_t& sat_idx,
 std::vector<size_t>* in_bins_idx)
{

  if(in_bins_idx)
  {
    if ( in_bins_idx->size() != m_SatBundles[sat_idx].GetBinPeriod().size() )
    {
      Z_LOG_ERROR(" The bin periods set size doesn't match the input size.") ;
      return;
    }
  }
  // else
  // {
  //   // if ( satPtrSet[sat_idx]->GetData().GetBins().size() != all_bin_periods.size() )
  //   // {
  //   //   Z_LOG_ERROR(" The bin periods set size doesn't match the dataset size.") ;
  //   //   return;
  //   // }
  // }


  // Cleaning up to avoid bugs
  // otherwise calling FitThreshold twice 
  // will keep growing the following vectors:
  if( thresh_limits_true.size()!=0 ||
      thresh_limits.size()!=0 ||  e_bins.size() != 0 )
  {
    Z_LOG_NOTE("Clearing the previous threshold fit results...") ;
    thresh_limits_true.clear() ;
    thresh_limits.clear() ;
    e_bins.clear() ;
  }
  //..........................................
  // ............ Creating a directory ............
  if (mkdir((m_SatBundles[sat_idx]->GetWrkDir()+"/Fit").Str().c_str(), ACCESSPERMS) == -1) 
  {
    Z_LOG_INFO("Directory '"+m_SatBundles[sat_idx]->GetWrkDir().Str()+"/Fit' wasn't created, because: "+strerror(errno)+".") ;  
  }
  else
    Z_LOG_INFO("Directory '"+m_SatBundles[sat_idx]->GetWrkDir().Str()+"/Fit' created."); 
  // .................................................
  if(in_bins_idx)
  {
    size_t j = 0 ;
    for(size_t i : *in_bins_idx)
    {
      // Bin b = satPtrSet[sat_idx]->GetData().GetBins()[i];

      // satPtrSet[sat_idx]->TimeBin(all_bin_periods[j++]);

      // size_t chops = satPtrSet[sat_idx]->GetSigShape(b.GetECenter().val).GetSize() -2  ;
      // // std::cout<<"\n Size - 2 = "<<chops<<" \n"<<std::flush;
      // b.divide(chops);
      // FitBin(b, b.GetECenter().val, sat_idx) ;
      FitBin(m_SatBundles[sat_idx].GetTimeBin(i,j++), sat_idx) ;
    }
  } 
  else
  {
    for(size_t i=0 ; i < m_SatBundles[sat_idx]->GetData().GetBins().size() ; ++i)
    {
      // Bin b = satPtrSet[sat_idx]->GetData().GetBins()[i];

      // satPtrSet[sat_idx]->TimeBin(all_bin_periods[i]);

      // size_t chops = satPtrSet[sat_idx]->GetSigShape(b.GetECenter().val).GetSize() -2  ;
      // // std::cout<<"\n Size - 2 = "<<chops<<" \n"<<std::flush;
      // b.divide(chops);
      // FitBin(b, b.GetECenter().val, sat_idx) ;
      FitBin(m_SatBundles[sat_idx].GetTimeBin(i,i), sat_idx) ;
    }
  }
  
  // Bin my_bin ;

}


//--------------------------------------------------------------
// Boosted Case
void Analysis::FitBoosted(Zaki::Math::Range<double> e_range)
{
  if(!modelPtr)
  {
    Z_LOG_ERROR("Model pointer hasn't been set yet, use 'SetModel'.") ;
    return ;
  }

  Zaki::String::Directory tmp_dir("") ;

  if(m_SatBundles.size() == 0 )
  {
    Z_LOG_ERROR("Satellite not attached!") ;
    return ;
  }

  if(m_SatBundles.size() == 1)
  {
    tmp_dir = m_SatBundles[0]->GetWrkDir() ;
  }
  else if (m_SatBundles.size() > 1)
  {  tmp_dir = wrk_dir + "/Combined" ;}
  else
  {
    Z_LOG_ERROR("Satellite not attached!") ;
    return ;
  }
  

  // ............ Creating a directory ............
  if (mkdir((tmp_dir + "/Fit").Str().c_str(), ACCESSPERMS) == -1) 
  {
    Z_LOG_INFO("Directory '"+tmp_dir.Str()+"/Fit' wasn't created, because: "+strerror(errno)+".") ;    
  }
  else
    Z_LOG_INFO("Directory '"+tmp_dir.Str()+"/Fit' created."); 
  // .................................................

  // Cleaning up to avoid bugs
  // otherwise calling FitThreshold twice 
  // will keep growing the following vectors:
  if( boost_fit_results.size()!=0 ||
      focus_bin_periods.size()!=0 ||  focus_bins.size() != 0 )
  {
    Z_LOG_NOTE("Clearing the previous boosted fit results...") ;
    focus_bins.clear() ;
    focus_bin_periods.clear() ;
    boost_fit_results.clear() ;
    boost_fit_val = 0 ;
  } 
  //....................................................

  // ROOT::Minuit2::MnMachinePrecision minuit_prec ;
  // minuit_prec.SetPrecision(1e-11) ;

  // Likelihood Function
  ROOT::Minuit2::LogLikeli  fcn ;

  // Fitting Parameters
  ROOT::Minuit2::MnUserParameters upar;

  // upar.SetPrecision(1e-8) ;

  // Signal Strength Parameter (mu)
  upar.Add("mu", 0.01, 1e-10);
  // upar.Add("mu", 1000, 0.00000001);

  upar.SetLowerLimit("mu", 0);

  //................................................................
  size_t focus_bins_size = 0 ;
  double mean_scale_factor = 0 ;
  double mean_exposure = 0 ;

  for(size_t sat_j=0 ; sat_j < m_SatBundles.size() ; ++sat_j)
  {
    focus_bins_size += m_SatBundles[sat_j]->GetData().GetBins().size() ;
    mean_scale_factor +=  m_SatBundles[sat_j]->GetNumScaling()/ m_SatBundles.size() ;
    mean_exposure +=  m_SatBundles[sat_j]->GetExpNorm()/ m_SatBundles.size() ;
  }

  // This Saves the allocation time
  focus_bins.reserve(focus_bins_size) ;
  focus_bin_periods.reserve(focus_bins_size) ;
  // _______________________1st Loop__________________________
  for(size_t sat_j=0 ; sat_j < m_SatBundles.size() ; ++sat_j)
  {
    // _______________________2nd Loop__________________________
    for(size_t i=0 ; i < m_SatBundles[sat_j]->GetData().GetBins().size() ; ++i)
    {
      Bin b = m_SatBundles[sat_j]->GetData().GetBins()[i];
      // ..............................
      // out of spectrum condition
      if( b.GetERange().min >= e_range.max || b.GetERange().max <= e_range.min )
        continue ; 
      
      focus_bins.push_back(b) ;
      focus_bin_periods.push_back(m_SatBundles[sat_j].GetBinPeriod()[i]) ;
      // ..............................

      // Time-binning
      m_SatBundles[sat_j]->TimeBin(m_SatBundles[sat_j].GetBinPeriod()[i]);

      // Getting the signal shape from Satellite 
      TH1F tmp_sig_hist = m_SatBundles[sat_j]->GetSigShape(b.GetECenter().val) ;

      // Dividing the bin into same pieces as the sigshape
      b.divide(tmp_sig_hist.GetSize() -2 );

      // Scaling the histogram by the energy-dependent 
      //  spectrum of the model
      tmp_sig_hist.Scale(modelPtr->GetSpectrum(b)) ;

      tmp_sig_hist.Scale(m_SatBundles[sat_j]->GetExpNorm()/mean_exposure) ;
      tmp_sig_hist.Scale(mean_scale_factor /m_SatBundles[sat_j]->GetNumScaling()) ;

      fcn.AddObsCounts(b) ;
      fcn.AddSigShape(tmp_sig_hist)  ;

      std::string bg_par_str = "bg_" + std::to_string(sat_j) 
                              + "_"+ std::to_string(i);  

      Zaki::Math::Quantity tmp_obs_q = b.GetTBinObsSet()[0] ;
      // double tmp_min_lim = ( tmp_obs_q.val > 2*tmp_obs_q.err > 0 ) ? 
      //                       tmp_obs_q.val - 2*tmp_obs_q.err : 0.5 ;

      // Adding the i-th background paramater with 
      // initial value 38 (could be anything reasonable)
      upar.Add(bg_par_str, tmp_obs_q.val, 1, 
              0, 4500);
      // upar.Add(bg_par_str, b.GetTBinObsSet()[0].val, 0.05,
      //          b.GetTBinObsSet()[0].val - b.GetTBinObsSet()[0].err,
      //          b.GetTBinObsSet()[0].val + b.GetTBinObsSet()[0].err);
      // upar.SetLowerLimit(bg_par_str, 0.2);
      // upar.SetUpperLimit(bg_par_str, 1000);
      // upar.SetLimits(bg_par_str, 1, 100 ) ;

    }
  // _________________________________________________________
  }
  // _________________________________________________________

  // Creat MIGRAD minimizer
  int strategy_counter = 2 ;
  ROOT::Minuit2::MnMigrad migrad(fcn, upar, strategy_counter);

  // migrad.SetPrecision(1e-8) ;

  // Minimize
  ROOT::Minuit2::FunctionMinimum min = migrad(4000);
  //.................
  // If minimization failed:
  while(!min.IsValid() && strategy_counter < 2) 
  {
    strategy_counter++ ;
    // try with higher strategy
    Z_LOG_ERROR("Migrad min is invalid, trying with strategy = "+ std::to_string(strategy_counter)+ ".");
    ROOT::Minuit2::MnMigrad migrad(fcn, upar, strategy_counter);
    min = migrad();
  }
  //.................

  // std::cout << "Minimum EDM = "<< min.Edm() << "\n" ;

  // output
  std::cout<<"minimum: "<<min<<std::endl;
  
  unsigned int Minos_max_f_calls = 4000 ;

  ROOT::Minuit2::MnMinos Minos(fcn, min);

  // MINOS errors
  // -2Log & chi2 have the same normalization:
  // Zaki::Math::NDimContLevel cont_level((int)upar.Params().size(), conf_level) ;
  // std::cout << "--> upar.Params().size(): " << (int)upar.Params().size() <<"\n" ;
  // fcn.SetErrorDef(cont_level.Solve().Up);
  fcn.SetErrorDef(ROOT::MathMore::chisquared_quantile(conf_level,
                  (int)upar.Params().size())) ;
  // fcn.SetErrorDef(pow(2.486,2));

  ROOT::Minuit2::MinosError e0 = Minos.Minos(0, Minos_max_f_calls); // get parameter 0's error
  // ROOT::Minuit2::MinosError e1 = Minos.Minos(1); // get parameter 1's error 

  boost_fit_results.reserve(focus_bins.size()+1) ;
  for(size_t i=0 ; i < focus_bins.size()+1 ; ++i)
  {
    boost_fit_results.push_back( Minos.Minos(i, Minos_max_f_calls)) ;
  }
  
  // boost_fit_val = (e0.Min() + e0.Upper() ) *
                  // (mean_scale_factor / mean_exposure) ;

  boost_fit_val = (e0.Min() + e0.Upper() ) * mean_scale_factor ;

  std::cout << "--> mean_scale_factor= " << mean_scale_factor << "\n";
  std::cout << "--> mean_exposure= " << mean_exposure << "\n";

  // ..............................................................
  // Creates Likelihood MnContours
  // for (size_t i = 0; i < focus_bins.size(); i++)
  // {
  //   ROOT::Minuit2::MnContours contour(fcn, min) ;

  //   fcn.SetErrorDef(ROOT::MathMore::chisquared_quantile(0.68, (int)upar.Params().size())) ;
  //   std::vector<std::pair<double, double>> pts_68 = contour(0,i+1, 150) ;

  //   fcn.SetErrorDef(ROOT::MathMore::chisquared_quantile(0.90, (int)upar.Params().size())) ;
  //   std::vector<std::pair<double, double>> pts_90 = contour(0,i+1, 150) ;
    
  //   fcn.SetErrorDef(ROOT::MathMore::chisquared_quantile(0.95, (int)upar.Params().size())) ;
  //   std::vector<std::pair<double, double>> pts_95 = contour(0,i+1, 150) ;
             
  //   TCanvas c("c", m_SatBundles[0]->GetName().c_str(), 1000, 800) ;
  //   c.SetGrid();                                                      
                             
  //   auto mg = new TMultiGraph();
    
  //   TGraph g68, g90, g95 ;     
    
  //   for (size_t i = 0; i < pts_68.size() ; i++)
  //   {
  //     g68.SetPoint(i, pts_68[i].first * m_SatBundles[0]->GetNumScaling() , pts_68[i].second) ;
  //   }                      
  //   g68.SetMarkerStyle(29) ;     
  //   g68.SetLineColor(kRed) ;
  //   g68.SetMarkerColor(kRed) ;

  //   for (size_t i = 0; i < pts_90.size() ; i++)
  //   {
  //     g90.SetPoint(i, pts_90[i].first * m_SatBundles[0]->GetNumScaling() , pts_90[i].second) ;
  //   }                      
  //   g90.SetMarkerStyle(29) ;     
  //   g90.SetLineColor(kGreen) ;
  //   g90.SetMarkerColor(kGreen) ;

  //   for (size_t i = 0; i < pts_95.size() ; i++)
  //   {
  //     g95.SetPoint(i, pts_95[i].first * m_SatBundles[0]->GetNumScaling() , pts_95[i].second) ;
  //   }                      
  //   g95.SetMarkerStyle(29) ;    
  //   g95.SetLineColor(kBlue) ;
  //   g95.SetMarkerColor(kBlue) ;

  //   mg->Add(&g68) ;
  //   mg->Add(&g90) ;
  //   mg->Add(&g95) ;
    
  //   mg->GetXaxis()->SetLimits(1e-5, 5e-3) ; 
  //   // mg->SetMinimum(0) ;
  //   // mg->SetMaximum(11) ;

  //   mg->GetXaxis()->SetMaxDigits(2) ;
  //   mg->Draw("ALP");      

  //   char out_file_char[150] ;
  //   sprintf(out_file_char, "M=%.0f/M_%.0f_B_%zu", modelPtr->GetDMMass(), modelPtr->GetDMMass(), i) ;
  //   c.SaveAs((wrk_dir+"/"+m_SatBundles[0]->GetName() + "/Fit/MnConts/"+ out_file_char  + ".pdf").Str().c_str()) ; 

  // }
  // // ..............................................................

  // // output
  // Z_LOG_INFO("95% Minos errors: ") ;
  // // std::cout<<e0<<"\n";

  // std::cout<< "\n * ------------------------------------------------------- * " ;
  // for (size_t i = 0; i < boost_fit_results.size() ; i++)
  // {
  //   std::cout<<boost_fit_results[i]<<"\n";
  // }
  // std::cout<< " * ------------------------------------------------------- * \n" ;


  // double tmp_sum = 0 ;
  // ------------------------------------------------------
  //    !!!           Temporary           !!!
  // ------------------------------------------------------
  // for (size_t i = 0; i < focus_bins.size() ; i++)
  // {
  //   Bin b2 = focus_bins[i] ;
  //   TH1F tmp_sig_hist = m_SatBundles[0]->GetSigShape(b2.GetECenter().val) ;

  //   b2.divide(tmp_sig_hist.GetSize() -2 );

  //   tmp_sig_hist.Scale(modelPtr->GetSpectrum(b2)) ;
  //   tmp_sig_hist.Scale(m_SatBundles[0]->GetExpNorm()/mean_exposure) ;
  //   tmp_sig_hist.Scale(mean_scale_factor /m_SatBundles[0]->GetNumScaling()) ;

  //   ROOT::Minuit2::LogLikeli  fcn2 ;
  //   fcn2.AddObsCounts(b2) ;
  //   fcn2.AddSigShape(tmp_sig_hist)  ;
  //   double tmp_val = fcn2({boost_fit_results[0].Min() + boost_fit_results[0].Upper(),
  //                                        boost_fit_results[i+1].Min() + boost_fit_results[i+1].Lower()}) ;
  //   tmp_val -= fcn2({boost_fit_results[0].Min(),
  //                                        boost_fit_results[i+1].Min()}) ;
  //   tmp_sum += tmp_val  ;
  //   std::cout << "\n" << i+1 << " ) " <<  tmp_val
  //                                        << "\n";
  // }
  // std::cout<< " Total =  "<< tmp_sum <<"\n" ;
  // std::cout<< " * ------------------------------------------------------- * \n" ;
  // // ------------------------------------------------------


  char tmp_char[150] ;
  sprintf(tmp_char, "\n -Lower Limit (95%%): %.2e\n -Upper Limit (95%%): %.2e",
          e0.Min() + e0.Lower(), e0.Min() + e0.Upper()) ;
  Z_LOG_INFO(tmp_char) ;

  sprintf(tmp_char, " ---> Minos error is invalid, # of calls for mu: %u", e0.NFcn()) ;
  if(!e0.IsValid()) Z_LOG_ERROR(tmp_char) ;
}

//-------------------------------------------------------------------------------
double Analysis::GetBoostLimit() const
{
  return boost_fit_val ;
}

//-------------------------------------------------------------------------------
void Analysis::SetConfidenceLevel(const double in_cl)
{
  conf_level = in_cl ;
}

//-------------------------------------------------------------------------------
// Generating the histograms in FitBin method for the threshold fit
double Analysis::Histgen(const size_t& sat_idx,
                         ROOT::Minuit2::MinosError e0,
                         ROOT::Minuit2::MinosError e1,
                         Bin b, TH1F sig_shape)
{

  double mu = e0.Min() + e0.Upper() ;
  double bg = e1.Min() ;
  size_t bin_num        = b.GetTBinChops()  ;
  int t_min             = 0                 ;   
  int t_max             = bin_num           ;

  TH1F o = b.GetTBinHist() ;

  TH1F tmp_sig_hist("sig", "Signal", bin_num, t_min, t_max);


  // Filling the signal histogram
  // Note the use of GetSize() and the fact that indexes of TH1F 
  // have 2 extra elements, at their first, and last position.
  for(int i=0; i<sig_shape.GetSize(); i++)
  {
    // Technically for i = 0 & i = sig_shape.GetSize() - 1
    // the element is '0'. 
    tmp_sig_hist.SetBinContent(i, sig_shape[i]*mu + bg)  ;  
  }

  TCanvas c("c", m_SatBundles[sat_idx]->GetName().c_str(), 1000, 600) ;
  c.SetGrid();
  // c.Divide(2) ;  
  // c.cd(1) ;
  // gPad->SetGrid() ;
  gStyle->SetOptStat(0);

  // Signal
  tmp_sig_hist.SetMarkerColor(kRed);
  tmp_sig_hist.SetMarkerStyle(5);
  tmp_sig_hist.SetLineColor(kRed);
  
  Zaki::Math::Range<double> tmp_obs = {b.GetTBinObsSet()[0].val - b.GetTBinObsSet()[0].err,
                     b.GetTBinObsSet()[0].val + b.GetTBinObsSet()[0].err } ;

  Zaki::Math::Range<double> tmp_sigbg = {tmp_sig_hist.GetMinimum(), tmp_sig_hist.GetMaximum()} ;

  Zaki::Math::Range<double> tmp_y = { tmp_sigbg.min< tmp_obs.min? tmp_sigbg.min*0.9 : tmp_obs.min*0.9 , 
                  tmp_sigbg.max > tmp_obs.max ? tmp_sigbg.max*1.1 : tmp_obs.max*1.1 };

  tmp_sig_hist.SetMaximum(tmp_y.max);
  tmp_sig_hist.SetMinimum(tmp_y.min);
  tmp_sig_hist.SetXTitle("Cycle") ;
  tmp_sig_hist.SetYTitle("Observed Counts") ;
  tmp_sig_hist.SetTitle("Constant Observed Events") ;

  tmp_sig_hist.Draw("LP") ;

  // Observed
  o.SetLineWidth(2) ;
  o.Draw("E1 same") ;

  TPaveText pl(t_max*0.7, tmp_y.max*0.8, t_max*0.98, tmp_y.max*0.98) ;

  char tmp_label[150] ;    
  sprintf(tmp_label, "E = [ %.1f - %.1f ] GeV", b.GetERange().min, b.GetERange().max) ;
  pl.AddText(tmp_label);

  sprintf(tmp_label, "#mu = %.3e #plus %.2e #minus %.2e", e0.Min(), e0.Upper(), abs(e0.Lower()) ) ;
  pl.AddText(tmp_label);

  sprintf(tmp_label, "bg = %.2f #plus %.2f #minus %.2f", e1.Min(), e1.Upper(), abs(e1.Lower()) ) ;
  pl.AddText(tmp_label);

  sprintf(tmp_label, "obs = %.2f #pm %.2f", b.GetTBinObsSet()[0].val, b.GetTBinObsSet()[0].err ) ;
  pl.AddText(tmp_label);

  pl.Draw();

  // Output file
  char out_file_char[150] ;    
  sprintf(out_file_char, "Fit/Thresh_%.0f_GeV.pdf", b.GetECenter().val) ;

  c.SaveAs((wrk_dir+"/"+m_SatBundles[sat_idx]->GetName() + "/"+ out_file_char).Str().c_str()) ;

  double muSig = 0 ;

  // Finding the overall limit
  for(int i=0; i<sig_shape.GetSize(); i++)
  {
    muSig += sig_shape[i]*mu ;
  }
  
  char tmp_char[100] ;
  sprintf(tmp_char, "--->>> 95%% C.L. on mu.S = %f", muSig) ;
  Z_LOG_INFO(tmp_char) ;

  return muSig ;
}

//-------------------------------------------------------------------------------
void Analysis::PlotThresh(const size_t& sat_idx, 
                          const Zaki::String::Directory& f_name)
{
  TCanvas c("c", m_SatBundles[sat_idx]->GetName().c_str(), 800, 600) ;
  c.SetGrid();
  
  gStyle->SetOptStat(0);

  TGraph gr(thresh_limits.size(), &e_bins[0], &thresh_limits[0]);
  // gr.SetMaximum(700);
  // gr.SetMinimum(50);
  gr.SetTitle("95% limit as a function of energy bin;E (GeV);#mu_{95}") ;
  gr.SetLineColor(kPink) ;
  gr.SetLineStyle(2) ;
  gr.SetMarkerColor(kPink) ;
  gr.SetMarkerSize(2) ;
  gr.SetMarkerStyle(29) ;
  gr.Draw("ALP");

  c.SaveAs((wrk_dir + "/" + m_SatBundles[sat_idx]->GetName()+ "/Fit/"+ f_name + ".pdf").Str().c_str()) ;
}

//--------------------------------------------------------------
void Analysis::PlotBoost(const Zaki::String::Directory& f_name)
{

  if(boost_fit_results.size() == 0)
  {
    Z_LOG_ERROR("Plotting failed, run FitBoosted first.") ;
    return ;
  }

  std::string tmp_name = "Combined" ;
  if(m_SatBundles.size() == 1)
    tmp_name = m_SatBundles[0]->GetName() ;

  TCanvas c("c", tmp_name.c_str(), 1000, 800) ;
  c.SetGrid();

  TGraphAsymmErrors gr;

  double max_y = 0 ;
  for(size_t i = 0 ; i<boost_fit_results.size()-1 ; ++i)
  {
    double e_h = focus_bins[i].GetERange().max -  focus_bins[i].GetECenter().val ;
    double e_l = focus_bins[i].GetECenter().val - focus_bins[i].GetERange().min;

    gr.SetPoint(i, focus_bins[i].GetECenter().val, boost_fit_results[i+1].Min()) ;
    gr.SetPointError(i, e_l, e_h, abs(boost_fit_results[i+1].Lower()),
                     boost_fit_results[i+1].Upper()) ;

    // For the location of plot legend
    double tmp_y =  boost_fit_results[i+1].Min() + boost_fit_results[i+1].Upper();
  
    if (max_y < tmp_y)
      max_y      = tmp_y;

  }

  gr.GetXaxis()->SetTitle("Characteristic Energy (#tilde{E})   [GeV]");
  gr.GetYaxis()->SetTitle("Background Fit Parameter");
  gr.SetTitle("Fit Results (Boosted)");
  gr.SetMarkerColor(4);
  gr.SetMarkerStyle(21);
  gr.Draw("ALP");

  // Plot legend
  double max_energy =focus_bins[focus_bins.size() - 1].GetERange().max;
  TPaveText pl(max_energy*0.2,
            max_y*0.88, 
            max_energy*0.98,
            max_y*0.98) ;  

  double mu_sc = 0 ;
  for(auto& satB : m_SatBundles)
    mu_sc += satB->GetNumScaling() / m_SatBundles.size(); 

  char tmp_label[150] ;  
  sprintf(tmp_label, "#mu = %.2e #plus %.2e #minus %.2e", 
          boost_fit_results[0].Min()*mu_sc,
          boost_fit_results[0].Upper()*mu_sc,
          abs(boost_fit_results[0].Lower())*mu_sc) ;
  pl.AddText(tmp_label);

  pl.Draw();

  c.SetLogx();

  // Output file  
  c.SaveAs((wrk_dir + "/" + tmp_name+ "/Fit"+ f_name + ".pdf").Str().c_str()) ;
}

//--------------------------------------------------------------
void Analysis::SetBinPeriod(const size_t& sat_idx,
                            const std::vector<int>& per)
{
  // all_bin_periods = per;
  m_SatBundles[sat_idx].SetBinPeriod(per) ;
  // for(int per : m_SatBundles[sat_idx].GetBinPeriod())
  //   std::cout << " Bin Period: " << per << "\n" ;
}

//--------------------------------------------------------------
std::vector<int> Analysis::GetBinPeriod(const size_t& sat_idx) const
{
  // return all_bin_periods;
  return m_SatBundles[sat_idx].GetBinPeriod() ;
}

//--------------------------------------------------------------
// Setting the model
void Analysis::SetModel(std::shared_ptr<Model> model_in)
{
  Z_LOG_INFO("Setting the input model...") ;
  modelPtr = model_in ;
  
  modelPtr->SetECut(GetECut()) ;
}

//--------------------------------------------------------------
/// Sets the plot mode
void Analysis::SetPlotMode(const PlotMode& in_mode) 
{
  plot_mode = in_mode ;
}

//--------------------------------------------------------------
/// Updates the input Member Function Contour Wrapper (MFCW)
void Analysis::UpdateMFCW(CONFIND::MemFuncContWrapper<Analysis, 
                          double (Analysis::*) (double, double)>* in_mfcw )
{
  switch (plot_mode)
  {
    case PlotMode::Ldec_Gann:
      in_mfcw->UpdateMemFunc(*this, &Analysis::ContFuncThresh_LG) ;
      break;
    
    case PlotMode::Mdm_Gann:
      in_mfcw->UpdateMemFunc(*this, &Analysis::ContFuncThresh_MG) ;
      break;

    case PlotMode::Ldec_Mdm:
      in_mfcw->UpdateMemFunc(*this, &Analysis::ContFuncThresh_LM) ;
      break;

    default:
      break;
  }
}

//--------------------------------------------------------------
/// Threshold analysis
/// If in_bins_idx is null it will fit all the satelite's bin
void Analysis::DoThreshMinuit(const size_t& sat_idx, 
  std::vector<size_t>* in_bins_idx, Zaki::Math::Grid2D* grid_in)
{
  PROFILE_FUNCTION() ;
  // Z_TIMER_SCOPE("Analysis::DoThresh") ;
  using namespace CONFIND ;

  if(GetSatellite(sat_idx)->GetData().size() == 0)
  {
    Z_LOG_ERROR("Satellite data is empty!") ;
    exit(EXIT_FAILURE) ;
  }
  
  //.........................................
  // Checking if the threshold fit has been done
  // if( thresh_limits_true.size() == 0)
  // {
    // Z_LOG_ERROR("Threshold fit hasn't been performed.") ;
  // Z_LOG_INFO("Performing the threshold fit...") ;
  FitThreshold(sat_idx, in_bins_idx) ;
  // }
  //.........................................

  // Initializing the model
  if(modelPtr)
    modelPtr->Init() ;
  else
  {
    Z_LOG_ERROR("Model is not set, use 'SetModel' first.") ;
    return;
  }
  
  // ...........................
  // Setting up the grid
  // ...........................
  // CONFIND::ContourFinder con    ;
  MemFuncContWrapper<Analysis, double (Analysis::*) (double, double)>* mfcwPtr;

  switch (plot_mode)
  {
    case PlotMode::Ldec_Gann:
      mfcwPtr = new MemFuncContWrapper<Analysis,
                double (Analysis::*) (double, double)>(*this, &Analysis::ContFuncThresh_LG) ;
      break;
    
    case PlotMode::Mdm_Gann:
      mfcwPtr = new MemFuncContWrapper<Analysis,
                double (Analysis::*) (double, double)>(*this, &Analysis::ContFuncThresh_MG) ;
      break;

    case PlotMode::Ldec_Mdm:
      mfcwPtr = new MemFuncContWrapper<Analysis, 
                double (Analysis::*) (double, double)>(*this, &Analysis::ContFuncThresh_LM) ;
      break;

    default:
      break;
  }

  // CONFIND::MemFuncContWrapper<Analysis, double (Analysis::*) (double, double)> mfcw(*this, &Analysis::ContFuncThresh_LG) ;
  // Zaki::Math::MemFuncWrapper<Analysis, double (Analysis::*) (double, double)> mfw(*this, &Analysis::ContFuncThresh) ;
  // con.SetMemFunc(&mfw);
  //..................
  // Checking the grid input
  if (grid_in)
  {
    Z_LOG_INFO("Using the grid input from 'Analysis'.") ;
    (*mfcwPtr)->SetGrid(*grid_in) ;
  }
  // Checking the grid in the model
  else if(modelPtr->GetTGrid(plot_mode))
  {
    Z_LOG_INFO("Using the grid from 'Model'.") ;
    (*mfcwPtr)->SetGrid(modelPtr->GetTGrid(plot_mode)->grid) ;
  }
  // Emitting errors and stopping the process
  else
  {
    Z_LOG_ERROR("Analysis aborted: No grid was input, and there is also none in the model!") ;
    return ;
  }
  //..................
  std::string tmp_name = "Combined" ; 
  if(m_SatBundles.size() == 1)
    tmp_name = m_SatBundles[sat_idx]->GetName() ;

  (*mfcwPtr)->SetWrkDir(wrk_dir + "/" + tmp_name + "/Fit") ;

  (*mfcwPtr)->SetThreads(8) ;

  // ...........................
  // Saving time by checking if the grid values are the
  //  same for various contours
  if(modelPtr->GetTGrid(plot_mode)->fixed)
  {
    std::vector<std::string> labels ;
    for (size_t i = 0; i < e_bins.size(); ++i)
    { 
      if(in_bins_idx)
        labels.emplace_back("T_"+std::to_string(GetBinPeriod(sat_idx)[i])) ;
      else
        labels.emplace_back("M_" +std::to_string((int)e_bins[i]*modelPtr->GetDecayProd().size())) ; 
    }

    (*mfcwPtr)->SetContVal(GetThreshLimits(), labels) ;
    (*mfcwPtr)->SetGridVals(CONFIND::ContourFinder::Mode::Fast) ;
  }
  else
  {
    for (size_t i = 0; i < thresh_limits_true.size(); ++i)
    {
      switch (plot_mode)
      {
        case PlotMode::Ldec_Gann:
          // Changing the DM mass in the model
          modelPtr->SetDMMass(e_bins[i]*modelPtr->GetDecayProd().size()) ;
        break;

        case PlotMode::Mdm_Gann:

        break;

        case PlotMode::Ldec_Mdm:

        break;

        default:
        break;
      }

      // // Recalculating other model parameters with the new DM mass
      // modelPtr->Init() ;

      if(in_bins_idx)
      {
        modelPtr->SetActiveBin(m_SatBundles[sat_idx]->GetData().GetBins()[(*in_bins_idx)[i]]) ;
        (*mfcwPtr)->SetContVal({thresh_limits_true[i]}, {"T_"+std::to_string(GetBinPeriod(sat_idx)[i])}) ;
      }
      else
      {
        modelPtr->SetActiveBin(m_SatBundles[sat_idx]->GetData().GetBins()[i]) ;
        (*mfcwPtr)->SetContVal({thresh_limits_true[i]}, {"M_"+std::to_string((int)e_bins[i]*modelPtr->GetDecayProd().size())}) ;
      }

      // Need to update mfw after any change
      // mfcwPtr->UpdateMemFunc(*this, &Analysis::ContFuncThresh) ;
      UpdateMFCW(mfcwPtr) ;

      // modelPtr->SetModelPars(dphot_params) ;
      // (*mfcwPtr)->SetContVal({thresh_limits_true[i]}, {std::to_string((int)e_bins[i]*modelPtr->GetDecayProd().size())}) ;
      (*mfcwPtr)->SetGridVals(CONFIND::ContourFinder::Mode::Fast) ;
    }
  }
  // ...........................

  // ...........................
  // Generating results
  // ...........................
  // (*mfcwPtr)->Print() ;
  // ............ Creating a directory ............
  if (mkdir((wrk_dir + "/" + tmp_name + "/Fit/Contour Points/Threshold").Str().c_str(), ACCESSPERMS) == -1) 
  {
    Z_LOG_INFO("Directory '"+ (wrk_dir + "/" + tmp_name + "/Fit/Contour Points/Threshold").Str()+"' wasn't created, because: "+strerror(errno)+".") ;    
  }
  else
    Z_LOG_INFO("Directory '"+(wrk_dir + "/" + tmp_name + "/Fit/Contour Points/Threshold").Str()+"' created."); 
  // .................................................
  (*mfcwPtr)->ExportContour("Contour Points/Threshold/"+ modelPtr->GetName() + "_threshold_cont",  Zaki::File::FileMode::Write) ; 

  size_t tmp_bins_size = 0 ; 
  double tmp_leg_x     = 0 ;
  std::string tmp_leg_label = "M_{DM} GeV";
  if(in_bins_idx) 
  { 
    tmp_bins_size = in_bins_idx->size() ;
    tmp_leg_x = 0.05 ;
    tmp_leg_label = "T (Day)" ;
  } 
  else            
  { 
    tmp_bins_size = m_SatBundles[sat_idx]->GetData().GetBins().size() ; 
  }
  // ...........................
  // Making the plot legend
  (*mfcwPtr)->MakeLegend(true, tmp_leg_label.c_str(), "user") ;
  (*mfcwPtr)->GetLegend()->SetX1(0.75 + tmp_leg_x) ; (*mfcwPtr)->GetLegend()->SetY1(0.10) ;
  (*mfcwPtr)->GetLegend()->SetX2(0.90) ; (*mfcwPtr)->GetLegend()->SetY2(0.15 + tmp_bins_size*0.03) ;
  (*mfcwPtr)->GetLegend()->SetTextSize(0.025) ;

  std::vector<std::string> tmp_labels ;
  char tmp[200] ;

  for (size_t i = 0; i < tmp_bins_size; i++)
  {
    if(in_bins_idx)
      sprintf(tmp, " %d", m_SatBundles[sat_idx].GetBinPeriod()[i]) ;
    else
      sprintf(tmp, "%.1e", e_bins[i]*modelPtr->GetDecayProd().size()) ;

    tmp_labels.emplace_back(tmp) ;
  }
    
  
  (*mfcwPtr)->SetLegendLabels(tmp_labels) ;
  // ...........................

  if(modelPtr->GetTGrid(plot_mode)->option.joined)
    (*mfcwPtr)->SetPlotConnected() ;
  
  (*mfcwPtr)->SetWidth(modelPtr->GetTGrid(plot_mode)->option.w) ;
  (*mfcwPtr)->SetHeight(modelPtr->GetTGrid(plot_mode)->option.h) ;

  sprintf(tmp,"%s - Threshold (%s)", modelPtr->GetName().c_str(),
          m_SatBundles[sat_idx]->GetName().c_str()) ;
  (*mfcwPtr)->Plot(modelPtr->GetName() + "_Threshold", tmp, 
           modelPtr->GetTGrid(plot_mode)->option.xLabel.c_str(), 
           modelPtr->GetTGrid(plot_mode)->option.yLabel.c_str()) ;
  // ...........................

  delete mfcwPtr;
}

//--------------------------------------------------------------
// Boosted analysis
void Analysis::DoBoost(Zaki::Math::Grid2D* grid_in)
{
  PROFILE_FUNCTION() ;
  // Z_TIMER_SCOPE("Analysis::DoBoost") ;

  using namespace CONFIND ;

  if(!modelPtr)
  {
    Z_LOG_ERROR("Model is not set, use 'SetModel' first.") ;
    exit(EXIT_FAILURE) ;
  }

  for(auto& satB : m_SatBundles)
  {
    if(satB->GetData().size() == 0)
    {
      Z_LOG_ERROR("Satellite: '" + satB->GetName() +"' has no data!") ;
      exit(EXIT_FAILURE) ;
    }
  }

  // Plot legend labels
  std::vector<std::string> tmp_labels ;

  // ...........................
  // Setting up the grid
  // ...........................
  // CONFIND::ContourFinder con    ;
  MemFuncContWrapper<Analysis, double (Analysis::*) (double, double)>* mfcwPtr;

  switch (plot_mode)
  {
    case PlotMode::Ldec_Gann:
      mfcwPtr = new MemFuncContWrapper<Analysis,
                double (Analysis::*) (double, double)>(*this, &Analysis::ContFuncBoost_LG) ;
      break;
    
    case PlotMode::Mdm_Gann:
      mfcwPtr = new MemFuncContWrapper<Analysis,
                double (Analysis::*) (double, double)>(*this, &Analysis::ContFuncBoost_MG) ;
      break;

    case PlotMode::Ldec_Mdm:
      mfcwPtr = new MemFuncContWrapper<Analysis, 
                double (Analysis::*) (double, double)>(*this, &Analysis::ContFuncBoost_LM) ;
      break;

    default:
      break;
  }
  // CONFIND::MemFuncContWrapper<Analysis, double (Analysis::*) (double, double)> mfcw(*this, &Analysis::ContFuncBoost) ;
  // Zaki::Math::MemFuncWrapper<Analysis, double (Analysis::*) (double, double)> mfw(*this, &Analysis::ContFuncBoost) ;
  // mfcw->SetMemFunc(&mfw);
  //..................
  // Checking the grid input
  if (grid_in)
  {
    Z_LOG_INFO("Using the grid input from 'Analysis'.") ;
    (*mfcwPtr)->SetGrid(*grid_in) ;
  }
  // Checking the grid in the model
  else if(modelPtr->GetBGrid(plot_mode))
  {
    Z_LOG_INFO("Using the grid from 'Model'.") ;
    (*mfcwPtr)->SetGrid(modelPtr->GetBGrid(plot_mode)->grid) ;
  }
  // Emitting errors and stopping the process
  else
  {
    Z_LOG_ERROR("No grid was input, and there is also none in the model!") ;
    Z_LOG_ERROR("Analysis aborted!") ;
    return ;
  }
  //..................
  std::string tmp_name = "Combined" ;
  if(m_SatBundles.size() == 1)
    tmp_name = m_SatBundles[0]->GetName() ;

  (*mfcwPtr)->SetWrkDir(wrk_dir + "/" + tmp_name + "/Fit") ;

  // Initializing the model
  modelPtr->Init() ;

  char tmp[200] ;
  for (size_t i = 0; i < modelPtr->GetDMMassSet().size(); i++)
  {

    // Changing the DM mass
    modelPtr->SetDMMass(modelPtr->GetDMMassSet()[i]) ;

    FitBoosted({GetECut(), modelPtr->GetDMMass()}) ;

    sprintf(tmp, "Boosted_Fit_%.0f_%.0f", GetECut(), modelPtr->GetDMMass()) ;
    PlotBoost(tmp);

    (*mfcwPtr)->SetContVal({GetBoostLimit()}, {std::to_string((int)modelPtr->GetDMMass())}) ;

    // Need to update mfw after any change
    // mfcw.UpdateMemFunc(*this, &Analysis::ContFuncBoost) ;
    UpdateMFCW(mfcwPtr) ;

    // ...........................
    (*mfcwPtr)->SetThreads(8) ;
    (*mfcwPtr)->SetGridVals(CONFIND::ContourFinder::Mode::Fast) ;

    sprintf(tmp, "%.1e", modelPtr->GetDMMass()) ;
    tmp_labels.emplace_back(tmp) ;

    // ...........................
    ResetContainers() ;
  }

  // ...........................
  // Making the plot legend
  (*mfcwPtr)->MakeLegend(true, "M_{DM}  [GeV]", "user") ;
  (*mfcwPtr)->GetLegend()->SetX1(0.75) ; (*mfcwPtr)->GetLegend()->SetY1(0.10) ;
  (*mfcwPtr)->GetLegend()->SetX2(0.90) ; (*mfcwPtr)->GetLegend()->SetY2(0.25) ;
  (*mfcwPtr)->GetLegend()->SetTextSize(0.025) ;
  (*mfcwPtr)->SetLegendLabels(tmp_labels) ;
  // ...........................
  if(modelPtr->GetBGrid(plot_mode)->option.joined)
    (*mfcwPtr)->SetPlotConnected() ;

  (*mfcwPtr)->SetWidth(modelPtr->GetBGrid(plot_mode)->option.w) ;
  (*mfcwPtr)->SetHeight(modelPtr->GetBGrid(plot_mode)->option.h) ;

  // ...........................
  // Generating results
  // ...........................
  // (*mfcwPtr)->Print() ;
  (*mfcwPtr)->ExportContour("Contour Points/Boosted/"+ modelPtr->GetName() + "_boosted_cont",  Zaki::File::FileMode::Write) ; 



  sprintf(tmp,"%s - Boosted (%s)", 
          modelPtr->GetName().c_str(), tmp_name.c_str()) ;  
  (*mfcwPtr)->Plot( modelPtr->GetName() + "_Boosted", 
            tmp, modelPtr->GetBGrid(plot_mode)->option.xLabel.c_str(), 
           modelPtr->GetBGrid(plot_mode)->option.yLabel.c_str()) ;

  delete mfcwPtr;
}

//--------------------------------------------------------------
void Analysis::ResetContainers()
{
  focus_bins.clear() ;
  focus_bin_periods.clear() ;
  boost_fit_results.clear() ;
}

//--------------------------------------------------------------
// Function for plotting contours
// Not to be confused with the same method in 'Model'.
// Boosted case, in (L_dec, G_ann) plane
double Analysis::ContFuncBoost_LG(double dec_len, double gamma)
{
  double out = modelPtr->ContFuncBoost_LG(dec_len, gamma) ;
  // out       *= satPtrSet[sat_idx]->GetKappa(modelPtr->GetDecayProd()) ;

  return out ;
}

//--------------------------------------------------------------
// Function for plotting contours
// Not to be confused with the same method in 'Model'.
// Boosted case, in (M_dm, G_ann) plane
// ! Not implemented yet !
double Analysis::ContFuncBoost_MG(double m_dm, double gamma)
{
  // double out = modelPtr->ContFuncBoost_MG(m_dm, gamma) ;
  // out       *= GetSatellite(0)->GetMultiplicity(modelPtr->GetDecayProd()) ;

  // return out ;
  return 1;
}

//--------------------------------------------------------------
// Function for plotting contours
// Not to be confused with the same method in 'Model'.
// Boosted case, in (L_dec, M_dm) plane
// ! Not implemented yet !
double Analysis::ContFuncBoost_LM(double dec_len, double m_dm)
{
  // double out = modelPtr->ContFuncBoost_LM(dec_len, m_dm) ;
  // out       *= GetSatellite(0)->GetMultiplicity(modelPtr->GetDecayProd()) ;

  // return out ;
  return 1;
}

//--------------------------------------------------------------
/// Function for plotting contours
/// Not to be confused with the same method in 'Model'.
/// Threshold case, in (L_dec, G_ann) plane
double Analysis::ContFuncThresh_LG(double dec_len, double gamma)
{
  double out = modelPtr->ContFuncThresh_LG(dec_len, gamma) ;
  out       *= GetSatellite(0)->GetMultiplicity(modelPtr->GetDecayProd()) ;

  return out ;
}

//--------------------------------------------------------------
/// Function for plotting contours
/// Not to be confused with the same method in 'Model'.
/// Threshold case, (M_dm, G_ann) plane
double Analysis::ContFuncThresh_MG(double m_dm, double gamma)
{
  double out = modelPtr->ContFuncThresh_MG(m_dm, gamma) ;
  out       *= GetSatellite(0)->GetMultiplicity(modelPtr->GetDecayProd()) ;

  return out ;
}

//--------------------------------------------------------------
/// Function for plotting contours
/// Not to be confused with the same method in 'Model'.
/// Threshold case, (L_dec, M_dm) plane
double Analysis::ContFuncThresh_LM(double L_dec, double M_dm)
{
  double out = modelPtr->ContFuncThresh_LM(L_dec, M_dm) ;
  out       *= GetSatellite(0)->GetMultiplicity(modelPtr->GetDecayProd()) ;

  return out ;
}

//==============================================================
