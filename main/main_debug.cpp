/* 

  This code is for for Finding the limit on signal strength 
  assuming no modulation in the background and observed events.

  - Last updated by Zaki Mar 13 2020

*/

#include <Zaki/File/VecSaver.hpp>

// Local headers
#include "DMSS/AMS.hpp"
#include "DMSS/Analysis.hpp"
#include "DMSS/Data.hpp"
#include "DMSS/GenericModel.hpp"
#include "DMSS/DarkPhoton.hpp"

#include "omp.h"
#include <sys/stat.h>

// ........................
// ZContour
class Z 
{
  public:
    Z() { std::cout << "Z: created( " << this << " )\n";}
    ~Z() {std::cout << "Z: destroyed( " << this << " )\n";}
    Z(const Z& other_z) 
    { std::cout << "Z: cpy_constructor from " << &other_z << " -> " << this << "\n";  
      Z_f = other_z.Z_f->Clone() ;
    }

    std::unique_ptr<Zaki::Math::Func2D> Z_f = nullptr ;
    void SetF2Dptr(Zaki::Math::Func2D* in_f) { Z_f = std::unique_ptr<Zaki::Math::Func2D>(in_f) ; }
    void SetPtr2Cell() ;
} ;

// ........................
// Cell
class C 
{
  public:
  C() { std::cout << "C: created( " << this << " )\n"; }
  ~C() { std::cout << "C: destroyed( " << this << " )\n"; }
  void f() ;
  Z* C_zptr = NULL ;
} ;
// ........................
void  C::f()
{
  double tmp = C_zptr->Z_f->Eval(1, 2) ; 
  std::cout << "C: C_tmp: "<< tmp << "\n" ;
} ;
// ........................
void Z::SetPtr2Cell() 
{
  // Z tmp_z[2] = {*this, *this} ;
  omp_set_num_threads(1) ;
  #pragma omp parallel
  {
    std::cout << "Z: SetPtr2Cell() : 56 " << "\n" ;
    Z tmp_z(*this) ;
    std::cout << "Z: SetPtr2Cell() : 58 " << "\n" ;
  C c ;  c.C_zptr = &tmp_z ; c.f() ;
  std::cout << "Z: SetPtr2Cell() : 60 " << "\n" ;
  }
} 
// ........................
class A
{
  public:
    A() { std::cout << "A: created (" << this << ")\n";}
    ~A() { std::cout << "A: destroyed (" << this << ")\n"; }
    A(const A& other) : a_var(other.a_var), z_obj(other.z_obj)
    { std::cout << "A: cpy_constructor from " << &other << " -> " << this << "\n"; }

    int a_var = 0 ;

    double a_f(double x, double y) { a_var++ ; std::cout<< "A: a_var: " << a_var <<"\n"; return x+y ;}

    Z z_obj ;
    void a_setMemFunc() 
    {
      Zaki::Math::MemFuncWrapper<A, double (A::*) (double, double)> w(*this, &A::a_f) ;
      z_obj.SetF2Dptr(&w) ;
      z_obj.SetPtr2Cell() ;
    }
};
// ........................

//************//
//    MAIN
//************//
int main() {

  std::string dir(__FILE__) ;
  dir = dir.substr(0, dir.find_last_of("/")) ;
  dir = dir.substr(0, dir.find_last_of("/")) ;

  // { Year, month, day, {hour, minute, seconds} } UTC
  // DMSolarSignal::Date start_1 = {2013, 10, 20, {7, 35, 0}} ;
  Zaki::Physics::Date Jun21   = {2016, 6, 21, {0, 0, 0}} ;


  // printf("%.5f\n", d.UnixTDay()) ;
  // std::cout << "Date: " << d << "\n" ;

  AMS s;
  // s.SetTimeDuration({{2015, 5, 20, {7, 35, 0}}, {2016, 5, 20, {7, 35, 0}}}) ;
  // std::cout << s.GetSunPos({2015, 10, 20, {7, 35, 0}}) << "\n" ;
  s.SetWrkDir(dir + "/results/debug") ; 


  Zaki::Physics::Date start = s.GetTimeDuration().start ;
  Zaki::Physics::Date end   = s.GetTimeDuration().end ;

  // s.SetYPR({0, 0, 60}) ;
  // s.SetOrbitInclination(0) ;
  s.Print() ;
  // std::cout << "GetSatPos: "<< s.GetSatPos(VERNAL_2016.UnixTDay()) << "\n" ;
  // std::cout << "GetSatPosEclip: "<< s.GetSatPos(VERNAL_2016.UnixTDay()).ToEcliptic(VERNAL_2016) << "\n" ;
  std::cout << "GetSunPos: "<< s.GetSunPos(Jun21.UnixTDay()) << "\n" ;
  std::cout << "GetSunPosEclip: "<< s.GetSunPos(Jun21.UnixTDay()).ToEcliptic(Jun21) << "\n" ;
  // std::cout << "GetNormGEI: "<< s.GetNormGEI(start.UnixTDay()) << "\n" ;
  // std::cout << "GetNormLocal: "<< s.GetNormLocal(start.UnixTDay()) << "\n" ;

  return 0 ;

  // s.PlotOrbit(start, end_1, "Orbit") ;

  s.PlotSunProj(start, end, "Sunshine_Projection") ;
  s.PlotSolarBeta(start, end, "Solar_Beta_Angle") ;

  s.EvaluateExposure(start, end) ;
  s.ExportExposure("AMS_Debug_Exp.dat", Zaki::File::FileMode::Write) ;
  s.PlotExposure("Exposure", 900) ;

  return 0;

  // Data d ;
  // d.SetName("CALET") ;
  // d.Input("data/CALET_e+e-.dat") ;
  // d.SetStates({DMSolarSignal::Electron, DMSolarSignal::Positron}) ;

  // double ExpTime = 24*3600*780 ; // seconds
  // double A_geom = 1040*1e-4 ; // in  units of m^2*sr

  // for (size_t i = 0; i < d.size(); i++)
  // {
  //   int count =  round((d[i].GetERange().max - d[i].GetERange().min)*ExpTime*A_geom*d[i].GetFlux().val)  ;
  //   double count_rel_err = sqrt( pow(d[i].GetFlux().sysErr, 2) + pow(d[i].GetFlux().statErr, 2)) / d[i].GetFlux().val ;

  //   std::cout << count << "\t" << round(sqrt(count)*10)/10.0 << "\t" << count*count_rel_err << "\n" ;
  // }

  // return 0 ;


  // A a ;
  // a.a_setMemFunc() ;
  
  // return 0 ;

  // test t1(1);
  // MemFuncWrapper<test, double (test::*) (double, double)> w1(t1, &test::f) ;
  // std::cout<< "\n &w1.fObj:  " << &w1.fObj << ",  &t1: " << &t1 << 
  // " ,  w1.fObj->var: " << w1.fObj->var << "\n" ;
  // MemFuncWrapper<test, double (test::*) (double, double)> a[3] = {w1, w1, w1} ;

  // a[0].fObj->var = 5 ; a[1].fObj->var = 6 ; a[2].fObj->var = 7 ;
  // for (size_t i = 0; i < 3; i++)
  // {
  //   std::cout<< &(a[i].fObj) << ", " << a[i].fObj->var << "\n" ;
  // }
  // std::cout<< "\n w1.fObj->var: " << w1.fObj->var << "\n" ;


  
  //.......................
  // Analysis
  //.......................
  Analysis a1 ;
  a1.SetName("Debug Analysis") ;
  a1.SetWrkDir(dir + "/results") ;

  //.......................
  // Satellite Input
  //.......................
  AMS ams     ;
  a1.AttachSatellite(&ams) ;

  // ams.Print() ;
  
  // ams.PlotSunProj({2012, 10, 20, {7, 35, 0}}, {2016, 10, 20, {7, 35, 0}}, "Sunshine_Projection") ;
  // ams.PlotBeta({2012, 10, 20, {7, 35, 0}}, {2016, 10, 20, {7, 35, 0}}, "Solar_Beta_Angle") ; 

  // ams.EvaluateExposure(0, 5*360);
  // ams.ExportExposure("data/AMS_1800_days_exp.dat", Zaki::File::FileMode::Write) ;
  
  ams.ImportExposure("data/AMS_1800_days_exp.dat") ;
  
  ams.NormalizeExposure() ;
  ams.PlotExposure("Exposure", 900) ;

  ams.TimeBin(30) ;
  ams.PlotSigShape("AMS_Sig_Shape") ;

  ams.GetData().Plot("AMS_Results") ;

  a1.SetBinPeriod({30, 30, 30, 45, 45, 90, 90, 180, 180, 360, 360, 360, 720}) ;

  // Boosted case
  a1.FitBoosted({50, 1000});
  a1.PlotBoost("AMS/Boosted_Fit");

  // Threshold case (bin-by-bin): Chi Squared Method
  a1.FitThreshold() ;
  a1.PlotThresh("Limits");
  
  // Threshold case (bin-by-bin): LogLikelihood Method
  // a1.FitThreshold("AMS/Log_Like") ;
  // a1.PlotThresh("AMS/Log_Like/Limits") ;

  
  // a1.DoGenericThresh() ;

  // a1.DoDarkPhotonThresh() ;

  return 0;
}
