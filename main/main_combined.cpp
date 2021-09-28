/* 

  This code is for for Finding the limit on signal strength 
  assuming no modulation in the background and observed events.

  - Last updated by Zaki June 30, 2020

*/

#include <Zaki/File/VecSaver.hpp>
#include <Zaki/Util/ObjObserver.hpp>

// Local headers
#include "DMSS/AMS.hpp"
#include "DMSS/DAMPE.hpp"
#include "DMSS/CALET.hpp"
#include "DMSS/Analysis.hpp"
#include "DMSS/GenericModel.hpp"
#include "DMSS/DarkPhoton.hpp"
#include "DMSS/Prog.hpp"

#define DO_PROFILING 0

//************//
//    MAIN
//************//
int main() {
  
  using namespace Zaki::Util ;

  Zaki::String::Directory dir(__FILE__) ;
  ObjManager::SetFile(dir.ParentDir() +"/results/Combined/Objs_combined.txt") ;
  
#if DO_PROFILING
  Instrumentor::BeginSession("Main", dir +"/results/Combined/Combined_Profile.json");
#endif

  LogManager::SetLogLevels(LogLevel::Info, LogLevel::Verbose) ;
  LogManager::SetLogFile(dir.ParentDir() +"/results/Combined/log_Combined_analysis.txt") ;
  LogManager::SetBlackWhite(true) ;


  {
    Z_TIMER_SCOPE("combined_main") ;

    //.......................
    // Analysis
    //.......................
    Analysis a1 ;
        
    a1.SetName("Combined Analysis") ;
    a1.SetWrkDir(dir.ParentDir() +"/results") ;
    //.......................
    // Satellite Input
    //.......................
    std::shared_ptr<AMS> ams = std::make_shared<AMS>() ;
    // std::shared_ptr<AMS> ams_2 = std::make_shared<AMS>() ;

   std::shared_ptr<CALET> calet = std::make_shared<CALET>() ;
//    std::shared_ptr<DAMPE> dampe = std::make_shared<DAMPE>() ;

    a1.AttachSatellite(ams) ;
    // ams_2->SetName("AMS_2") ;
    // a1.AttachSatellite(ams_2) ;
   a1.AttachSatellite(calet) ;
//    a1.AttachSatellite(dampe) ;
    
    // ....................................................
    // Exposure
    ams->ImportExposure("AMS_Exposure_May2011_Nov2017.dat") ;
    ams->NormalizeExposure() ;
    ams->TimeBin(79) ;

    // ams_2->ImportExposure("AMS_Exposure_May2011_Nov2017.dat") ;
    // ams_2->NormalizeExposure() ;
    // ams_2->TimeBin(79) ;
    
   calet->ImportExposure("CALET_Exp_Oct2015_Nov2017.dat") ;
   calet->NormalizeExposure() ;
   calet->TimeBin(30) ;
//
//    dampe->ImportExposure("DAMPE_Exposure_Dec2015_Jun2017.dat") ;
//    dampe->NormalizeExposure() ;
//    dampe->TimeBin(10) ;

    a1.SetBinPeriod(0, {10, 10, 10, 15, 15, 15, 30, 30, 30, 30, 79, 79, 2*79}) ;
    // a1.SetBinPeriod(1, {10, 10, 10, 15, 15, 15, 30, 30, 30, 30, 79, 79, 2*79}) ;
   a1.SetBinPeriod(1, {20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 39, 39, 39, 39, 39, 130, 156, 195, 390, 390}) ;
//    a1.SetBinPeriod(2, {53, 53, 53, 53, 53, 53, 106, 106, 106, 265, 265, 265, 265, 530, 530, 530, 530}) ;
    
    a1.SetECut(50) ;
    
    // Boosted case
    // a1.FitBoosted({50, 1000}) ; // E_Cut = 50 GeV, M_DM = 1000 GeV
    // a1.PlotBoost("Boosted_Fit");

    // Threshold case (bin-by-bin): Chi Squared Method
//    a1.FitThreshold() ;
    // a1.PlotThresh("Limits") ;


    //........ Generic Model .............
   std::shared_ptr<GenericModel> genMod = std::make_shared<GenericModel>() ;
   a1.SetModel(genMod) ;

//    std::vector<size_t> idx_set (a1.GetBinPeriod().size(), 0) ;
    
    // Threshold Analysis:
//     a1.DoThresh(0) ;

    // Boosted Analysis:
     genMod->SetDMMass({500, 1000, 1500, 2000}) ;
     a1.DoBoost()  ;

    //.....................................

    //........ Dark Photon ................
    // std::shared_ptr<DarkPhoton> darkMod = std::make_shared<DarkPhoton>() ;
    // a1.SetModel(darkMod) ;

    // Threshold Analysis:
  //  a1.DoThresh(0) ;

    // Boosted Analysis:
    // darkMod->SetDMMass({500, 1000, 1500, 2000}) ;
//    darkMod->SetDMMass({1000}) ;
    // a1.DoBoost()  ;

    //.....................................
  }
  
#if DO_PROFILING
  Instrumentor::EndSession();
#endif

  return 0;
}
