/* 
*/

#include <Zaki/String/Directory.hpp>

// #include 

// Local headers
#include "DMSS/Fitter.hpp"
#include "DMSS/AMS.hpp"


//************//
//    MAIN
//************//
int main() {

  Zaki::String::Directory dir(__FILE__) ;

  std::shared_ptr<AMS> ams = std::make_shared<AMS>() ;


  // Getting the start & end times of operation
  auto start = ams->GetTimeDuration().start ;
  auto end   = ams->GetTimeDuration().end   ;

  // Setting the time stamp for the file names
  char tmp[150] ;
  sprintf(tmp, "_%s_%s", start.StrForm("#b#Y").c_str(),
          end.StrForm("#b#Y").c_str()) ;
  std::string time_stamp(tmp) ;

  ams->SetWrkDir(dir.ParentDir() +"/results") ;

  ams->ImportExposure("AMS_Exposure" + time_stamp + ".dat") ;

  Fitter fitter ;
  fitter.SetWrkDir(dir.ParentDir() +"/results") ;


  std::vector<double> time_bins = {10, 10, 10, 15, 15, 15, 30, 30, 30, 30, 79, 79, 2*79} ;
  std::vector<double> old_results = {12.52+0.08, 11.72+0.05, 11.11+0.02, 12.47+0.03, 11.40+0.06, 10.22+0.08, 42.09+0.09, 36.27+0.09,
                                     30.25+0.09, 24.66+0.09, 15.48+0.09, 9.42+0.06, 3.68+0.09 } ;
  
  for(size_t i=0 ; i < ams->GetData().size() ; i++)
  {
    Bin b = ams->GetData().GetBins()[i] ;

    ams->TimeBin(time_bins[i]);

    size_t chops = ams->GetSigShape(b.GetECenter().val).GetSize() -2  ;
    b.divide(chops);

    TH1F sig_shape = ams->GetSigShape(b.GetECenter().val) ;

    fitter.AddObsCounts(b) ;
    fitter.AddSigShape(sig_shape)  ;
  }

  std::vector<double> results = fitter.Fit() ;

  for(size_t i=0 ; i < results.size() ; i++)
  {
    std::cout << "\nScaled Result [" << i << "]: "<< results[i]* ams->GetNumScaling() / ams->GetExpNorm() << "\t" ;
    std::cout << "Old Result= " << (ams->GetNumScaling() / ams->GetExpNorm())*old_results[i] << "\n" ;
  }
  
//  std::cout << "ams->GetNumScaling() / ams->GetExpNorm()= " << ams->GetNumScaling() / ams->GetExpNorm() << "\n" ;
  // std::cout << "mu_old * ams->GetNumScaling() / ams->GetExpNorm()= " << 12.60*ams->GetNumScaling() / ams->GetExpNorm() << "\n" ;

  return 0;
}
