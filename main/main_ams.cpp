/* 

  This code is for for Finding the limit on signal strength 
  assuming no modulation in the background and observed events.

  - Last updated by Zaki Aug 17, 2020

*/

#include <Zaki/File/VecSaver.hpp>
#include <Zaki/Util/ObjObserver.hpp>

// Local headers
#include "DMSS/AMS.hpp"
#include "DMSS/Analysis.hpp"
#include "DMSS/GenericModel.hpp"
#include "DMSS/DarkPhoton.hpp"
#include "DMSS/Prog.hpp"


#define EVAL_EXP 0
#define DO_PROFILING 0
#define FUTURE 1

//************//
//    MAIN
//************//
int main() {
  
  using namespace Zaki::Util ;

  std::string folder_name = "" ;
#if FUTURE
  folder_name = "_2028";
#endif

  Zaki::String::Directory dir(__FILE__) ;
  
#if DO_PROFILING
  Instrumentor::BeginSession("Main", dir +"/results/AMS"+folder_name+"/AMS_Profile.json");
#endif

  ObjManager::SetFile(dir.ParentDir() +"/results/AMS"+folder_name+"/Objs_ams.txt") ;
  LogManager::SetLogLevels(LogLevel::Info, LogLevel::Verbose) ;
  LogManager::SetLogFile(dir.ParentDir() +"/results/AMS"+folder_name+"/log_AMS_analysis.txt") ;
  // LogManager::SetBlackWhite(true) ;

    Z_TIMER_SCOPE("ams_main") ;

    //.......................
    // Analysis
    //.......................
    Analysis a1 ;
        
    a1.SetName("AMS Analysis") ;
    a1.SetWrkDir(dir.ParentDir() +"/results") ;
    //.......................
    // Satellite Input
    //.......................
    std::shared_ptr<AMS> ams = std::make_shared<AMS>() ;


#if FUTURE
    ams->SetName("AMS_2028");
    // ams->SetTimeDuration({{2011, 5, 19, {0, 0, 0}}, {2021, 3, 27, {0, 0, 0}}}) ; 
    ams->SetTimeDuration({{2011, 5, 19, {0, 0, 0}}, {2027, 10, 22, {0, 0, 0}}}) ; 
#else
    ams->SetTimeDuration(ams->GetActualTimeDuration()) ; 
#endif

    ams->Print() ;  // Prints the Satellite info

    a1.AttachSatellite(ams) ;

    // std::cout << "\n "<< ams->GetTimeDuration().duration() <<", Ratio=" <<
    //  ams->GetTimeDuration().duration() / ams->GetActualTimeDuration().duration() 
    //  << "\n" ;
    // return 0 ;

    //  Getting the start & end times of operation
    auto start = ams->GetTimeDuration().start ;
    //     {2017, 11, 13, {0, 0, 0}}
    auto end   = ams->GetTimeDuration().end ;
    // Zaki::Physics::Date end = {2020, 11, 13, {0, 0, 0}}  ;

    // Comparing with DarkSunshine paper:
//     Zaki::Physics::Date start(2012, 1, 1, {0, 0, 0}) ;
//     Zaki::Physics::Date end(2012, 12, 31, {0, 0, 0}) ;


    // Setting the time stamp for the file names
    char tmp[150] ;
    sprintf(tmp, "_%s_%s", start.StrForm("#b#Y").c_str(),
            end.StrForm("#b#Y").c_str()) ;
    std::string time_stamp(tmp) ;

    // ......................................................................
    // Checking reference points from Ref:
    // "CALCULATION OF THE SUN, MOON AND ISS POSITIONS"
    // NICOLA TOMASSETTI, INFN
    // DMSolarSignal::Date ref_date1 = {2011, 7, 15, {0, 0 ,0}} ;
    // DMSolarSignal::Date ref_date2 = {2012, 1, 15, {0, 0 ,0}} ;
    // std::cout << "\n RAAN(ref1): " <<  ams->GetRAAN(ref_date1.UnixTDay()) << "\n" ;
    // std::cout << "\n RAAN(ref2): " <<  ams->GetRAAN(ref_date2.UnixTDay()) << "\n" ;
    // ......................................................................

    ams->PlotSunProj(start, end,  "Sunshine_Projection" + time_stamp) ;
    ams->PlotSolarBeta(start, end, "Solar_Beta_Angle" + time_stamp) ;

    // ....................................................
    // Exposure
#if EVAL_EXP
    ams->EvaluateExposure(start, end);
    ams->ExportExposure("AMS_Exposure" + time_stamp + ".dat", Zaki::File::FileMode::Write) ;
    return 0 ;
#else
    ams->ImportExposure("AMS_Exposure" + time_stamp + ".dat") ;
#endif
// return 0 ;
    ams->NormalizeExposure() ;
    ams->PlotExposure("Exposure" + time_stamp, 900) ;
    // ....................................................

    ams->TimeBin(79) ;
    ams->PlotSigShape("AMS_Sig_Shape" + time_stamp ) ;

    ams->GetData().Plot("AMS_Results" + time_stamp) ;
// return 0 ;


#if FUTURE
   a1.SetBinPeriod(0, {10, 10, 10, 10, 15, 15, 15, 15, 15, 30, 45, 80, 3*80}) ;
//  a1.SetBinPeriod(0, {10, 10, 10, 10, 10, 10, 10, 13, 20, 30, 45, 80, 3*80}) ;

  // Time bins for comparison
  // a1.SetBinPeriod(0, {10, 25, 50, 100, 150, 1185}) ;
#else
  // 13 bins
  // Factors of 2370 : { 2 × 3 × 5 × 79 }
  // a1.SetBinPeriod(0, {10, 10, 10, 10, 15, 15, 15, 15, 15, 30, 79, 79, 3*79}) ;
  
  // Time bins for all the bins 
  a1.SetBinPeriod(0, {10, 10, 10, 10, 15, 15, 15, 15, 15, 30, 45, 80, 3*80}) ;
  // a1.SetBinPeriod(0, {10, 10, 10, 10, 10, 10, 10, 13, 20, 30, 45, 80, 3*80}) ;

  // a1.SetBinPeriod(0, {1185, 1185, 1185, 1185, 1185, 1185, 1185, 1185, 1185, 1185, 1185, 1185, 1185}) ;

  // Time bins for comparison
  // a1.SetBinPeriod(0, {10, 25, 50, 100, 150, 1185}) ;
#endif

  a1.SetECut(50) ;
  
  // Boosted case
  // a1.FitBoosted({50, 1000}) ; // E_Cut = 50 GeV, M_DM = 1000 GeV
  // a1.PlotBoost("Boosted_Fit");

  // Threshold case (bin-by-bin)
  //  a1.FitThreshold() ;
  // a1.PlotThresh("Limits") ;

  //........ Generic Model .............
  std::shared_ptr<GenericModel> genMod = std::make_shared<GenericModel>() ;
  a1.SetModel(genMod) ;
  
  // std::vector<size_t>  idx_set = {0, 0, 0, 0, 0, 0} ;
  // Threshold Analysis:
  // a1.DoThresh(0, &idx_set) ;
  // a1.DoThreshMinuit(0, &idx_set) ;
  // a1.PlotParamSpaceThresh(0, false) ;

  // Boosted Analysis:
  // genMod->SetDMMass({200, 300, 400, 500, 750, 1000, 1500, 2000, 3000, 4000, 5000}) ;
  // genMod->SetDMMass({200, 225, 250, 275, 300, 325, 350, 375, 400, 425, 450, 475, 
                      // 500, 550, 600, 650, 700, 
                      // 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200, 1250,
                      // 1300, 1350,  1400, 1500,
                      // 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 
                      // 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400,
                      // 3500, 3600, 3700, 3800, 3900, 4000, 4100, 4200, 4300, 4400,
                      //  4500, 4600, 4700, 4800, 4900, 5000}) ;

  genMod->SetDMMass({  325, 350, 425, 600, 700, 800 }) ;

  // a1.PlotParamSpaceBoost(0, false) ;
  // return 0 ;

  a1.DoBoost()  ;
  //.....................................

  // return 0 ;

  //........ Dark Photon ................
  std::shared_ptr<DarkPhoton> darkMod = std::make_shared<DarkPhoton>() ;
  a1.SetModel(darkMod) ;

  // a1.PlotParamSpace() ; 

  // Threshold Analysis:
  // a1.DoThresh(0) ;

  // Boosted Analysis:
  // darkMod->SetDMMass({200, 300, 400, 500, 750, 1000, 1500, 2000, 3000, 4000, 5000}) ;
  darkMod->SetDMMass({  325, 350, 425, 600, 700, 800 }) ;

  a1.DoBoost()  ;
return 0 ;

  //.....................................
  
#if DO_PROFILING
  Instrumentor::EndSession();
#endif

  return 0;
}
