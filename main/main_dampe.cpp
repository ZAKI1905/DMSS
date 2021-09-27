/* 

  This code is for for Finding the limit on signal strength 
  assuming no modulation in the background and observed events.

  - Last updated by Zaki Mar June 4, 2020

*/

#include <Zaki/File/VecSaver.hpp>
#include <Zaki/Util/ObjObserver.hpp>

// Local headers
#include "DMSS/DAMPE.hpp"
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
  ObjManager::SetFile(dir.ParentDir() +"/results/DAMPE"+folder_name+"/Objs_dampe.txt") ;

#if DO_PROFILING
  Instrumentor::BeginSession("DAMPE",  dir +"/results/DAMPE"+folder_name+"/DAMPE_Profile.json");        // Begin session 
#endif

  LogManager::SetLogLevels(LogLevel::Info, LogLevel::Verbose) ;
  LogManager::SetLogFile(dir.ParentDir() +"/results/DAMPE"+folder_name+"/log_DAMPE_analysis.txt") ;
  // LogManager::SetBlackWhite(true) ;

  //.......................
  // Analysis
  //.......................
  Analysis a1;
  a1.SetName("DAMPE Analysis");
  a1.SetWrkDir(dir.ParentDir() + "/results") ;
  
  //.......................
  // Satellite Input
  //.......................
  std::shared_ptr<DAMPE> dampe = std::make_shared<DAMPE>() ;

  #if FUTURE
    dampe->SetName("DAMPE_2020");
    // dampe->SetTimeDuration({{2011, 5, 19, {0, 0, 0}}, {2021, 3, 27, {0, 0, 0}}}) ; 
    dampe->SetTimeDuration({{2015, 12, 27, {0, 0, 0}}, {2020, 11, 30, {0, 0, 0}}}) ; 
#else
    dampe->SetTimeDuration(dampe->GetActualTimeDuration()) ; 
#endif

  a1.AttachSatellite(dampe) ;
  dampe->Print();
  
  // Getting the start & end times of operation 
  auto start = dampe->GetTimeDuration().start ;
  auto end   = dampe->GetTimeDuration().end   ; 

  // Setting the time stamp for the file names
  char tmp[150] ;
  sprintf(tmp, "_%s_%s", start.StrForm("#b#Y").c_str(),
          end.StrForm("#b#Y").c_str()) ;
  std::string time_stamp(tmp) ;

  dampe->PlotSunProj(start, end,  "Sunshine_Projection" + time_stamp) ;
  dampe->PlotSolarBeta(start, end, "Solar_Beta_Angle" + time_stamp); 

  // ....................................................
  // Exposure 


#if EVAL_EXP
  dampe->EvaluateExposure(start, end);
  dampe->ExportExposure("DAMPE_Exposure" + time_stamp + ".dat", Zaki::File::FileMode::Write) ;
#else
  dampe->ImportExposure("DAMPE_Exposure" + time_stamp + ".dat") ;
#endif

  dampe->NormalizeExposure() ;
  dampe->PlotExposure("DAMPE_Exposure" + time_stamp, 900) ;
  // ....................................................

  dampe->TimeBin(10);

  dampe->PlotSigShape("DAMPE_Sig_Shape" + time_stamp);

  dampe->GetData().Plot("DAMPE_Results" + time_stamp);


#if FUTURE
  // a1.SetBinPeriod(0, {10, 10, 10, 10, 10, 10, 10, 10 ,10, 10, 10,
                      // 18, 18, 18, 18, 30, 30, 90, 100, 100, 225, 225,
                      // 450, 450, 900, 900, 900, 900}) ;

  a1.SetBinPeriod(0, {10, 10, 10, 10, 10, 10, 10, 10 ,10, 10, 10,
                      18, 18, 18, 18, 30, 30, 90, 100, 100, 265, 265,
                      265, 265, 530, 530, 530, 530}) ;
#else
  // 530 days
  // a1.SetBinPeriod(0, {53, 53, 53, 53, 53, 53, 106, 106, 106, 265, 265,
  //                     265, 265, 530, 530, 530, 530}) ;
  a1.SetBinPeriod(0, {10, 10, 10, 10, 10, 10, 10, 10 ,10, 10, 10,
                      18, 18, 18, 18, 30, 30, 90, 100, 100, 265, 265,
                      265, 265, 530, 530, 530, 530}) ;
#endif
  
  a1.SetECut(50) ; 

  // Boosted case
  // a1.FitBoosted();
  // a1.PlotBoost("DAMPE/Boosted_Fit");

  // Threshold case (bin-by-bin)
  // a1.FitThreshold() ;
  // a1.PlotThresh("Limits");
// return 0 ;
  //........ Generic Model .............
  std::shared_ptr<GenericModel> genMod = std::make_shared<GenericModel>() ;
  a1.SetModel(genMod) ;

  a1.PlotParamSpaceThresh(0, false) ; 
// return 0 ;

  // Threshold Analysis:
  // Zaki::Math::Grid2D thresh_dampe_grid = {{{5e4, 1e13}, 200, "Log"}, {{1e-6, 9e-1}, 200, "Log"}} ;
  // a1.DoThresh(0, nullptr, nullptr) ;

  // Boosted Analysis:
  // genMod->SetDMMass({200, 225, 250, 275, 300, 325, 350, 375, 400, 425, 450, 475, 
                    // 500, 550, 600, 650, 700, 
                    // 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200, 1250,
                    // 1300, 1350, 1400, 1500,
                    // 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 
                    // 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400,
                    // 3500, 3600, 3700, 3800, 3900, 4000, 4100, 4200, 4300, 4400,
                    //  4500, 4600, 4700, 4800, 4900, 5000}) ;
  // genMod->SetDMMass({4200}) ;

  // a1.PlotParamSpaceBoost(0, false) ; 
  // a1.DoBoost()  ;
  // return 0 ;
  //.....................................

  //........ Dark Photon .. .............
  std::shared_ptr<DarkPhoton> darkMod = std::make_shared<DarkPhoton>() ;
  a1.SetModel(darkMod) ;

  a1.PlotParamSpaceThresh(0, true) ; 
return 0 ;
  // Threshold Analysis:
  a1.DoThresh(0) ;

  // Boosted Analysis:
  // darkMod->SetDMMass({200, 300, 400, 500, 750, 1000, 1500, 2000, 3000, 4000, 5000}) ;
  // darkMod->SetDMMass({4200}) ;
  // a1.PlotParamSpaceBoost(0, true) ;
  // a1.DoBoost()  ;
  //.....................................

#if DO_PROFILING
  Instrumentor::EndSession();
#endif    

  return 0;

}
