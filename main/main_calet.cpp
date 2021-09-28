/* 

  This code is for for Finding the limit on signal strength 
  assuming no modulation in the background and observed events.

  - Last updated by Zaki Mar June 4, 2020

*/


#include <Zaki/File/VecSaver.hpp>
#include <Zaki/Util/ObjObserver.hpp>

// Local headers
#include "DMSS/CALET.hpp"
#include "DMSS/Analysis.hpp"
#include "DMSS/GenericModel.hpp"
#include "DMSS/DarkPhoton.hpp"

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
  folder_name = "_2020";
#endif

  Zaki::String::Directory dir(__FILE__) ;

  ObjManager::SetFile(dir.ParentDir() +"/results/CALET"+folder_name+"/Objs_calet.txt") ;

#if DO_PROFILING
  Instrumentor::BeginSession("CALET",  dir +"/results/CALET"+folder_name+"/CALET_Profile.json");        // Begin session 
#endif
  LogManager::SetLogLevels(LogLevel::Info, LogLevel::Verbose) ;
  LogManager::SetLogFile(dir.ParentDir() +"/results/CALET"+folder_name+"/log_CALET_analysis.txt") ;
  // LogManager::SetBlackWhite(true) ;

  //.......................
  // Analysis
  //.......................
  Analysis a1 ;
  a1.SetName("CALET Analysis") ;
  a1.SetWrkDir(dir.ParentDir() +"/results") ;

  //.......................
  // Satellite Input
  //.......................
  std::shared_ptr<CALET> calet = std::make_shared<CALET>() ;

#if FUTURE
    calet->SetName("CALET_2020");
    // ams->SetTimeDuration({{2011, 5, 19, {0, 0, 0}}, {2021, 3, 27, {0, 0, 0}}}) ; 
    calet->SetTimeDuration({{2015, 10, 13, {0, 0, 0}}, {2020, 9, 16, {0, 0, 0}}}) ; 
#else
    calet->SetTimeDuration(calet->GetActualTimeDuration()) ; 
#endif

  a1.AttachSatellite(calet) ;

  // Zaki::Physics::DateInterval inter = {{2015, 10, 13, {0, 0, 0}}, {2020, 9, 16, {0, 0, 0}}} ;
  // std::cout << "inter: " << inter.duration() << "\n" ;
  calet->Print() ;

  // Getting the start & end times of operation 
  auto start = calet->GetTimeDuration().start ;
  auto end   = calet->GetTimeDuration().end   ;

  // ......................................................................
  // Checking reference points from Ref:
  // "CALCULATION OF THE SUN, MOON AND ISS POSITIONS"
  // NICOLA TOMASSETTI, INFN 
  // DMSolarSignal::Date ref_date1 = {2011, 7, 15, {0, 0 ,0}} ;
  // DMSolarSignal::Date ref_date2 = {2012, 1, 15, {0, 0 ,0}} ;
  // std::cout << "\n RAAN(ref1 = 20.5609): " <<  calet->GetRAAN(ref_date1.UnixTDay()) << "\n" ;
  // std::cout << "\n RAAN(ref2 = 174.6112): " <<  calet->GetRAAN(ref_date2.UnixTDay()) << "\n" ;
  // ......................................................................

  // Setting the time stamp for the file names
  char tmp[150] ;
  sprintf(tmp, "_%s_%s", start.StrForm("#b#Y").c_str(),
          end.StrForm("#b#Y").c_str()) ; 
  std::string time_stamp(tmp) ;

  calet->PlotSunProj(start, end, "Sunshine_Projection" + time_stamp) ;
  calet->PlotSolarBeta(start, end, "Solar_Beta_Angle" + time_stamp) ; 

  // ....................................................
  // Exposure 

#if EVAL_EXP
  calet->EvaluateExposure(start, end);
  calet->ExportExposure("CALET_Exp" + time_stamp + ".dat", Zaki::File::FileMode::Write) ;
#else
  calet->ImportExposure("CALET_Exp" + time_stamp + ".dat") ;
#endif

  calet->NormalizeExposure() ;
  calet->PlotExposure("Exposure" + time_stamp, 900) ;
  // ....................................................

  calet->TimeBin(30) ;
  calet->PlotSigShape("CALET_Sig_Shape" + time_stamp) ;

  calet->GetData().Plot("CALET_Results" + time_stamp) ;


#if FUTURE
  a1.SetBinPeriod(0, {18, 18, 18, 18, 18, 18, 18, 18, 18, 20, 20, 20, 
                      30, 36, 36, 36, 60, 60, 90, 300, 360}) ;
#else
  // 21 bins 
  // Factors of 780: 2,3,4,5,6,10,12,13,15,20,26,30,39,52,60,65,78,130,156,195,260,390,780
  a1.SetBinPeriod(0, {18, 18, 18, 18, 18, 18, 18, 18, 18, 20, 20, 20, 
                      30, 36, 36, 36, 60, 60, 90, 300, 360}) ;
  // a1.SetBinPeriod(0, {20, 20, 20, 26, 26, 26, 26, 30, 30, 30, 39, 39, 39,
  //                     52, 52, 78, 130, 156, 195, 260, 390}) ;
#endif


  a1.SetECut(50) ; 

  // Boosted case
  // a1.FitBoosted({50, 1000}) ; // E_Cut = 50 GeV, M_DM = 1000 GeV
  // a1.PlotBoost("Boosted_Fit");

  // Threshold case (bin-by-bin): Chi Squared Method
  // a1.FitThreshold() ;
  // a1.PlotThresh("Limits") ;

  //........ Generic Model .............
  std::shared_ptr<GenericModel> genMod = std::make_shared<GenericModel>() ;
  a1.SetModel(genMod) ;

  // a1.PlotParamSpaceThresh(0, true) ; 
// return 0 ;

  // Threshold Analysis:
//  Zaki::Math::Grid2D thresh_calet_grid = {{{5e4, 1e13}, 200, "Log"}, {{3e-6, 9e-1}, 200, "Log"}} ;
//  a1.DoThresh(0, nullptr, nullptr) ;

  // Boosted Analysis:
  // genMod->SetDMMass({200, 225, 250, 275, 300, 325, 350, 375, 400, 425, 450, 475, 
                    // 500, 550, 600, 650, 700, 
                    // 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200, 1250,
                    // 1300, 1350, 1400, 1500,
                    // 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 
                    // 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400,
                    // 3500, 3600, 3700, 3800, 3900, 4000, 4100, 4200, 4300, 4400,
                    //  4500, 4600, 4700, 4800, 4900, 5000}) ;
  genMod->SetDMMass({ 4100, 4200, 4300, 4400, 4500, 4600, 4700, 4800, 4900 }) ;
  // a1.PlotParamSpaceBoost(0, true) ;
  // return 0 ;
  a1.DoBoost()  ;
  // return 0 ;
  //.....................................

  //........ Dark Photon .. .............
  std::shared_ptr<DarkPhoton> darkMod = std::make_shared<DarkPhoton>() ;
  a1.SetModel(darkMod) ;

  // Threshold Analysis:
  // a1.DoThresh(0) ;
  // a1.PlotParamSpaceThresh(0, true) ; 
  // return 0 ;

  // Boosted Analysis:
  // darkMod->SetDMMass({200, 300, 400, 500, 750, 1000, 1500, 2000, 3000, 4000, 5000}) ;
  darkMod->SetDMMass({ 4100, 4200, 4300, 4400, 4500, 4600, 4700, 4800, 4900 }) ;
  a1.DoBoost()  ;
  //.....................................

#if DO_PROFILING
  Instrumentor::EndSession();
#endif                    

  return 0;
}
