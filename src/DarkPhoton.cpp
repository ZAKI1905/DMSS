/*
  DarkPhoton class

*/


#include <gsl/gsl_integration.h>

// ROOT
#include <TF1.h>
#include <TCanvas.h>
#include <TStyle.h>
// #include <Math/WrappedFunction.h>
// #include <Math/Integrator.h>
#include <Math/Interpolator.h>

// Dependencies
#include <Zaki/Math/GSLFuncWrapper.hpp>
#include <Zaki/File/CSVIterator.hpp>
#include <Zaki/Vector/Vector_Basic.hpp>

// Local headers
#include "DMSS/DarkPhoton.hpp"

//==============================================================
using namespace Zaki::Physics ;
//--------------------------------------------------------------
// Constructor
DarkPhoton::DarkPhoton() 
{
  // logger.SetUnit("DarkPhoton") ;
  SetName("Dark Photon") ;
  SetElemSet() ; 

  size_t tmp_Res = 300; 

  //....................................................................//
  //                  Setting the boosted grids                         //
  //....................................................................//
  SetBGrid(PlotMode::Ldec_Gann , 
           {{{1.5e-3, 1e-1}, tmp_Res, "Log"},
            {{1e-10, 3e-6}, tmp_Res, "Log"}},
           {1300, 1000, "m_{A'} [GeV]", "#epsilon", false},
           false) ;
  SetBGrid(PlotMode::Ldec_Mdm , 
            {{{10*ELECTRON_M_GEV, 1}, tmp_Res, "Log"},
            {{100, 4000}, tmp_Res, "Linear"}},
           {1000, 1000, "m_{A'} [GeV]", "M_{DM} [GeV]", false},
           false) ;
  SetBGrid(PlotMode::Mdm_Gann , 
            {{{100, 4000}, tmp_Res, "Linear"},
            {{3e-10, 1e-6}, tmp_Res, "Log"}},
            {1000, 1000, "M_{DM} [GeV]", "#epsilon", false},
            false) ;
  //....................................................................//

  //....................................................................//
  //                  Setting the threshold grids                       //
  //....................................................................//
  SetTGrid(PlotMode::Ldec_Gann , 
           {{{2.001*ELECTRON_M_GEV, 4*ELECTRON_M_GEV}, tmp_Res, "Log"},
            {{3e-10, 1e-5}, tmp_Res, "Log"}},
           {1300, 1000, "m_{A'} [GeV]", "#epsilon", false},
           false) ;
  SetTGrid(PlotMode::Ldec_Mdm , 
            {{{2.001*ELECTRON_M_GEV, 4*ELECTRON_M_GEV}, tmp_Res, "Log"},
            {{100, 4000}, tmp_Res, "Log"}},
           {1000, 1000, "m_{A'} [GeV]", "M_{DM} [GeV]", false},
           false) ;
  SetTGrid(PlotMode::Mdm_Gann , 
            {{{100, 4000}, tmp_Res, "Log"},
            {{3e-10, 1e-5}, tmp_Res, "Log"}},
            {1000, 1000, "M_{DM} [GeV]", "#epsilon", false},
            false) ;
  //....................................................................//

  // // Setting the boosted grid
  // SetBGrid({{{10*ELECTRON_M_GEV, 1}, tmp_Res, "Log"},
  //         {{3e-10, 1e-6}, tmp_Res, "Log"}}) ;

  // // Setting the threshold grid
  // SetTGrid({{{2.01*ELECTRON_M_GEV, 
  //             10*ELECTRON_M_GEV},
  //           tmp_Res, "Log"}, {{3e-10, 6e-6}, tmp_Res, "Log"}}) ;
  
  // Setting the plot options:
  //  {width, height, x-Label, y-label, connected}
  // SetPlotOptions({1000, 1000, "m_{A'} [GeV]", "#epsilon", false}) ;

  // DPhotPar dphot_params ;
  // Excact value doesn't matter
  DPhotPar dphot_params ;
  dphot_params.eps = 1e-10 ;
  dphot_params.mDM = 2000  ;
  dphot_params.mDP = 0.1   ;
  SetModelPars(dphot_params) ;

  Z_LOG_NOTE(("DarkPhoton constructor called for: " + PtrStr()).c_str()) ;

}

//--------------------------------------------------------------
// Destructor
DarkPhoton::~DarkPhoton() 
{
  Z_LOG_NOTE(("DarkPhoton destructor called for: " + PtrStr()).c_str()) ;
}

//--------------------------------------------------------------
// Assignment operator 
// DarkPhoton& DarkPhoton::operator = (const DarkPhoton &t) 
// {
//   return *this;
// }

//--------------------------------------------------------------
// Copy constructor 
// DarkPhoton::DarkPhoton(const DarkPhoton &t)  
// {
//   label = t.label ;
//   set_label_flag = t.set_label_flag ;

//   dark_model = t.dark_model ;
//   set_model_par_flag = t.set_model_par_flag ;

//   DP_mass_br = t.DP_mass_br ;
//   set_mass_br_flag = t.set_mass_br_flag ;

//   element_set = t.element_set ;
//   set_elem_set_flag = t.set_elem_set_flag ;

//   fix_relic_flag = t.fix_relic_flag ;
//   num_decay_products = t.num_decay_products ;

//   br_data_file = t.br_data_file ;
//   // ROOT::Math::Interpolator inter ;
// }

//--------------------------------------------------------------
// IClone method
DarkPhoton* DarkPhoton::IClone() const
{
  DarkPhoton* tmp_ptr = new DarkPhoton(*this) ;

  std::stringstream ss;
  ss << static_cast<const void*>(tmp_ptr);

  Z_LOG_NOTE(("DarkPhoton::Clone() called from: " + PtrStr() + " to: " + ss.str()).c_str()) ;
  
  return tmp_ptr;
}

//--------------------------------------------------------------
// void DarkPhoton::SetLabel(const std::string& label_in) 
// {
//   set_label_flag = true ;
//   label = label_in ;
// }

//--------------------------------------------------------------
// std::string DarkPhoton::GetLabel() 
// {
//   if (set_label_flag)
//     return label ;
//   else 
//     return "[NOT SET]" ;
// }

// //--------------------------------------------------------------
// void DarkPhoton::SetDecayProd(const std::vector<DMSolarSignal::Particle>& in)
// {
//   decay_products = in ;
// }

// //--------------------------------------------------------------
// std::vector<DMSolarSignal::Particle> DarkPhoton::GetDecayProd() const
// {
//   return decay_products ;
// }

//--------------------------------------------------------------
// double DarkPhoton::GetDecayProdE() const
// {
//   if (!set_model_par_flag)
//   {
//     Z_LOG_ERROR("Model not set, use 'SetModelPars(DPhotPar x)' first.") ;
//     return -1 ;
//   }

//   return dark_model.mDM / decay_products.size() ;
// }

//--------------------------------------------------------------
void DarkPhoton::SetModelPars(DPhotPar& model_in)
{
  set_model_par_flag = true ;
  dark_model = model_in ;
}

//--------------------------------------------------------------
DPhotPar DarkPhoton::GetModelPars() const
{
  if (!set_model_par_flag)
    Z_LOG_ERROR("Model not set, use 'SetModelPars(DPhotPar x)' first.") ;

  return dark_model ;
}

//--------------------------------------------------------------
// Velocity of thermal DM in the sun
double DarkPhoton::GetSunDMVel() const
{
  if (!set_model_par_flag)
    Z_LOG_ERROR("Model not set, use 'SetModelPars(DPhotPar x)' first.") ;

  return sqrt(2.0*SUN_T_GEV / dark_model.mDM) ;
}
//--------------------------------------------------------------
// Approximate Sommerfeld formula (no resonances)
double DarkPhoton::GetApproxSomm() const
{
  if (!set_model_par_flag)
    Z_LOG_ERROR("Model not set, use 'SetModelPars(DPhotPar x)' first.") ;

  if( dark_model.alphaX() / GetSunDMVel() < 1e-10)
    return 1;

  double num = M_PI * dark_model.alphaX() / GetSunDMVel() ;
  double denom = 1-exp(-M_PI * dark_model.alphaX() / GetSunDMVel()) ;

  return num / denom ;
}

//--------------------------------------------------------------
// Sommerfeld (with resonances)  
// Ref: [ arXiv:1302.3898 eqs (34 -- 35) ]
void DarkPhoton::EvalExactSomm()
{
  if (!set_model_par_flag)
    Z_LOG_ERROR("Model not set, use 'SetModelPars(DPhotPar x)' first.") ;

  // if( dark_model.alphaX() / GetSunDMVel() < 1e-10)
  //   sommerfeld_factor = 1 ;

  // ROOT::Math::WrappedMemFunction<DarkPhoton, double ( DarkPhoton::* ) (double)>
  //  f(*this, &DarkPhoton::SommIntegrand);

  // ROOT::Math::Integrator ig(f) ;

  // ig.SetRelTolerance(0.01) ;
  double x_max_lim = 0.0005 ;
  // while ( f(x_max_lim) < 1e-01 ) x_max_lim -= 0.00001;

  // sommerfeld_factor = ig.Integral(0, x_max_lim) ;

  gsl_integration_workspace *w = gsl_integration_workspace_alloc(200);
  double err;

  // DarkPhoton local_model ;
  // local_model.Init() ;

  Zaki::Math::GSLFuncWrapper<DarkPhoton, double (DarkPhoton::*)(double)> 
    Fp(this, &DarkPhoton::SommIntegrand);     

  gsl_function F = static_cast<gsl_function> (Fp) ; 

  gsl_integration_qag(&F, 0, x_max_lim, 1e-6, 1e-6, 200, 1, w, &sommerfeld_factor, &err);
  gsl_integration_workspace_free(w);

  found_exact_somm_flag = true ;

}

//--------------------------------------------------------------
double DarkPhoton::GetExactSomm() const
{
  if (!found_exact_somm_flag)
  {
    Z_LOG_ERROR("Exact Sommerfeld factor hasn't been evaluated yet, use 'EvalExactSomm' first!") ;
    Z_LOG_WARNING(" ==> Returning approximate Sommerfeld instead.") ;
    return GetApproxSomm() ;
  }

  return sommerfeld_factor ;
}

//--------------------------------------------------------------
// Sommerfeld Integrand (with resonances)  
// Ref: [ arXiv:1302.3898 eqs (34 -- 35) ]
double DarkPhoton::SommIntegrand(double v)
{
  if ( v == 0 ) return 0 ;

  double v0 = GetSunDMVel() ;
  double a  = v / 2.0 / dark_model.alphaX() ;
  double c  = 6 * dark_model.alphaX() * dark_model.mDM / dark_model.mDP / pow(M_PI, 2) ;

  double Ss = M_PI / a;
  if ( 2*M_PI*a*c < 50 )
  {
    Ss       *= sinh(2*M_PI*a*c) ;
    Ss       *= 1.0 / ( cosh(2*M_PI*a*c) - cos(2*M_PI*sqrt(c - pow(a*c, 2))) ) ;
  }

  double Sav = 4*M_PI*pow(v, 2) ;
  Sav       *= 1.0 / pow(2*M_PI*pow(v0, 2), 3.0 / 2) ;
  Sav       *= exp(-0.5*pow(v / v0, 2)) * Ss;

  return Sav;
}

//--------------------------------------------------------------
// Sommerfeld Enhancement Factor 
double DarkPhoton::GetSommerfeld() const
{
  if ( exact_somm_flag )
    return GetExactSomm() ;
  else
    return GetApproxSomm() ;
}

//--------------------------------------------------------------
// DM annihilation cross section 
double DarkPhoton::GetAnnXSec() const
{
  if (!set_model_par_flag)
    Z_LOG_ERROR("Model not set, use 'SetModelPars(DPhotPar x)' first.") ;

  double fac_1 = pow( dark_model.gX, 4 ) / ( 16*M_PI*pow(dark_model.mDM, 2) ) ;

  double fac_2_num   = pow(1 - pow(dark_model.mDP/dark_model.mDM, 2) , 3.0/2) ;
  double fac_2_denom = pow(1 - 0.5*pow(dark_model.mDP/dark_model.mDM, 2),2) ;

  return fac_1 * fac_2_num / fac_2_denom ;
}

//--------------------------------------------------------------
// Fix gX assuming the correct relic abundance
void DarkPhoton::FixRelic() 
{
  if (!set_model_par_flag)
    Z_LOG_ERROR("Model not set, use 'SetModelPars(DPhotPar x)' first.") ;

  double qx_4 = SIG_V_REL*pow(CM_2_GEV, 3) / SEC_2_GEV ;

  qx_4 *= 16*M_PI*pow(dark_model.mDM, 2) ;

  qx_4 *= pow(1 - 0.5*pow(dark_model.mDP/dark_model.mDM, 2), 2) ;

  qx_4 *= 1.0 / pow(1 - pow(dark_model.mDP/dark_model.mDM, 2), 3.0/2) ;

  dark_model.gX = sqrt(sqrt(qx_4)) ;

  fix_relic_flag = true ;

  if ( exact_somm_flag )
    // #pragma omp critical
    { EvalExactSomm() ; }
}

//--------------------------------------------------------------
double DarkPhoton::GetCAnn() const
{
  double tmp   = dark_model.mDM*NEWTON_G_GEV*SUN_RHO_GEV ;
  tmp         *= 1.0 / (3.0* SUN_T_GEV) ;

  double c_ann = GetAnnXSec() ;
  c_ann       *= pow( tmp, 3.0/2);

  return c_ann;
}

//--------------------------------------------------------------
/*
  Setting the elements info taken from Table 1 of Ref:

  "A Poor Particle Physicist Cookbook for 
  Neutrinos from Dark Matter annihilations in the Sun"

  Link: https://arxiv.org/pdf/1312.6408.pdf 


  M.Cirelli, G.Corcella, A.Hektor, G.Hütsi, M.Kadastik, P.Panci, M.Raidal, F.Sala, A.Strumia,
  "PPPC 4 DM ID: A Poor Particle Physicist Cookbook for Dark Matter Indirect Detection'',
  arXiv 1012.4515, JCAP 1103 (2011) 051.
  Erratum: JCAP 1210 (2012) E01.

*/
void DarkPhoton::SetElemSet()
{
  set_elem_set_flag = true ;

  // Element e ;
  // e.sym  = "He" ;
  // e.Z    = 2 ;
  // e.A    = 4.002602 ; // in atomic unit
  // e.I    = 2*pow(10, 46) ;
  // e.v2   = 0.000017574529017920766 ;

  element_set.reserve(8) ;

  element_set.emplace_back("He", 2, 4.002602, 2*1e46, 0.000017574529017920766) ;

  element_set.emplace_back("N", 7, 14.0067, 2.84*1e43, 0.000017934976348385784) ;

  element_set.emplace_back("O", 8, 15.999, 7.34*1e43, 0.00001587280264793272) ;

  element_set.emplace_back("Ne", 10, 20.1797, 1.03*1e43, 0.000015908558331337145) ;

  element_set.emplace_back("Mg", 12, 24.305, 3.1*1e42, 0.000015819240309503544) ;

  element_set.emplace_back("Si", 14, 28.0855, 1.98*1e42, 0.000015819240309503544) ;

  element_set.emplace_back("S", 16, 32.065, 5.73*1e41, 0.00001587280264793272) ;

  element_set.emplace_back("Fe", 26, 55.845, 8.87*1e40, 0.000015819240309503544) ;
}

//--------------------------------------------------------------
// Elastic scattering for an element in the Sun
double DarkPhoton::GetElasticSig(Zaki::Physics::Element e) const
{
  double term = pow(e.Z, 2)*pow(e.mu(),2) ;
  term       *= dark_model.alphaX()*pow(dark_model.eps, 2) ;
  term       *= ALPHA_EM ;

  term *= 1.0 / pow(2*e.v2*pow(e.mu(), 2) + pow(dark_model.mDP, 2), 2) ;

  return 16 * M_PI * term;
}

//--------------------------------------------------------------
// DM capture rate for a specific element in the Sun
double DarkPhoton::GetCapRate(Zaki::Physics::Element e) const
{
  double cap_rate = 4*M_PI*F_Sun*DM_RHO_GEV / pow(dark_model.mDM, 2) ;
  
  cap_rate       *= GetElasticSig(e)*e.mu()*e.I;

  return cap_rate ;
}

//--------------------------------------------------------------
// DM capture rate from all elements in the Sun
double DarkPhoton::GetCapRateTot() const
{
  double tot_cap_rate = 0 ;

  for (size_t i = 0; i < element_set.size(); ++i)
  {
    tot_cap_rate += GetCapRate(element_set[i]) ;
  }

  return tot_cap_rate; 
}

//--------------------------------------------------------------
// Thermal equilibrium time in the Sun
double DarkPhoton::GetEqTau() const
{
  return 1.0 / sqrt(GetSommerfeld()*GetCAnn()*GetCapRateTot()); 
}

//--------------------------------------------------------------
// Annihilation rate (at equilibrium) in the Sun
double DarkPhoton::GetAnnRate() const
{
  return 0.5 * GetCapRateTot() * pow( tanh(SUN_AGE_GEV / GetEqTau()), 2); 
}

//--------------------------------------------------------------
double DarkPhoton::GetDecLenT() const
{
  if (dark_model.mDP <= 2*ELECTRON_M_GEV)
  {
    Z_LOG_WARNING("The dark photon decay is kinematically not possible.") ;
    return INFINITY ; // Some very large number
  }

  double dec_len = SUN_R_GEV*1; // Dark photon Br -> e is 1.

  dec_len       *= pow(1.1*1e-9/dark_model.eps, 2);
  dec_len       *= dark_model.mDM / dark_model.mDP / 3e6 ;
  dec_len       *= 1e-3 / dark_model.mDP ;
  dec_len       *= 2.32072e5 / sqrt(1 - pow(2*ELECTRON_M_GEV/dark_model.mDP, 2)) ;

  return dec_len;
}

//--------------------------------------------------------------
double DarkPhoton::GetDecLenB() const
{
  double dec_len = SUN_R_GEV*GetDPhotonBr();

  dec_len       *= pow(1.1*1e-9/dark_model.eps, 2);
  dec_len       *= dark_model.mDM / dark_model.mDP / 1e3 ;
  dec_len       *= 0.1 / dark_model.mDP ;

  return dec_len;
}

//--------------------------------------------------------------
double DarkPhoton::GetDecayProbT() const
{
  double p  = exp( - SUN_R_GEV / GetDecLenT() ) ;
  p        -= exp( - EARTH_2_SUN_GEV / GetDecLenT() ) ;
  
  return p;
}

//--------------------------------------------------------------
double DarkPhoton::GetDecayProbB() const
{
  double p  = exp( - SUN_R_GEV / GetDecLenB() ) ;
  p        -= exp( - EARTH_2_SUN_GEV / GetDecLenB() ) ;

  p        *= 1 - GetECut() / dark_model.mDM ;
  
  return p;
}

//--------------------------------------------------------------
// void DarkPhoton::SetECut(double e_cut_in)
// {
//   energy_cut = e_cut_in ;
// }

// //--------------------------------------------------------------
// double DarkPhoton::GetECut() const
// {
//   return energy_cut;
// }

//--------------------------------------------------------------
void DarkPhoton::ImportDarkPhotonBr(const Zaki::String::Directory& file_name)
{
  std::ifstream     file(file_name.Str());

  if (file.fail())
  {
    Z_LOG_ERROR(("Branching ratio file '" + file_name.Str() + "' cannot be opened!").c_str()) ;
    return;
  }
  else
  {
    Z_LOG_INFO(("Branching ratio file '" + file_name.Str() + "' imported.").c_str()) ;
  }

  // Reading the input file
  for(Zaki::File::CSVIterator loop(file); loop != Zaki::File::CSVIterator(); ++loop)
  {
    if (Zaki::Vector::Exists(static_cast<double>(std::stof((*loop)[0])), m_vals))
      continue ;
    DP_mass_br.push_back({std::stof((*loop)[0]), std::stof((*loop)[1]) }) ;
    m_vals.push_back(std::stof((*loop)[0])) ;
    br_vals.push_back(std::stof((*loop)[1])) ;
  }

  set_mass_br_flag = true ;
}

//--------------------------------------------------------------
double DarkPhoton::DarkPhotonBrWrapper(double *x, double *par) const
{
  if(*x < 2*ELECTRON_M_GEV)
    return 0;
  if(*x < 0.1810672687926081)
    return 1;
  if(*x > 1.6871421901704333)
    return 0.23903238192988588;
  
  // if (! set_mass_br_flag )
  // {
  //   std::cout<<"ERROR: Branching ratio data not imported!"<<"\n";
  //   std::cout<<"Importing Branching ratio from '"<<br_data_file<<"'...\n";
  //   ImportDarkPhotonBr(br_data_file) ;
  //   return -1 ;
  // }
  ROOT::Math::Interpolator inter ;
  inter.SetData(m_vals, br_vals);
  return inter.Eval(*x) ;
}

//--------------------------------------------------------------
double DarkPhoton::GetDPhotonBr() const
{
  if( dark_model.mDP < 2*ELECTRON_M_GEV )
    return 0 ;

  if(dark_model.mDP < 0.1810672687926081)
    return 1;

  if(dark_model.mDP > 1.6871421901704333)
    return 0.23903238192988588;

  // if (! set_mass_br_flag )
  // {
  //   std::cout<<"ERROR: Branching ratio data not imported!"<<"\n";
  //   std::cout<<"Importing Branching ratio from '"<<br_data_file<<"'...\n";
  //   ImportDarkPhotonBr() ;
  //   // return -1 ;
  // }

  ROOT::Math::Interpolator inter ;
  inter.SetData(m_vals, br_vals);
  return inter.Eval(dark_model.mDP) ;
}

//--------------------------------------------------------------
// Plot dark photon branching ratio
void DarkPhoton::DPhotonBrPlot(double m_1, double m_2, const std::string& f_name) const
{
  TCanvas c("c", "Dark Photon Branching Ratio", 1000, 600) ;
  c.SetGrid();
  gPad->SetGrid() ;
  gStyle->SetOptStat(0);

  TF1 BrPlot("BrPlot", this, &DarkPhoton::DarkPhotonBrWrapper, m_1, m_2, 0);
  BrPlot.SetTitle("Dark Photon Branching Ratio; m_{A'} [GeV]; Br") ;
  // BetaPlot.SetNpx(3500);

  
  // double max_y = BrPlot.GetMaximum(t_1, t_2) + 5 ;
  // double min_y = BrPlot.GetMinimum(t_1, t_2) - 5 ;

  // BrPlot.SetMaximum(max_y);
  // BrPlot.SetMinimum(min_y); 
  BrPlot.Draw() ;

  c.SaveAs((wrk_dir + "/" + f_name + ".pdf").Str().c_str()) ;
}

//--------------------------------------------------------------
// Function for plotting contours (inherited from base class)
// Threshold case in (m_A', eps) plane
double DarkPhoton::ContFuncThresh_LG(double x, double y)
{
  PROFILE_FUNCTION() ;
  dark_model.eps = y ;
  dark_model.mDP = x ;

  // Just for the flag mostly
  // SetModelPars(dark_model) ;

  FixRelic();

  double  model_factor = 1.0 ; // Br = 1
  model_factor        *= GetDecayProbT()     ;
  // model_factor        *= GetDecayProdNum()   ;
  
  // ( X X ==> fi fi ) Annihilation
  model_factor        *= 2*GetAnnRate()        ;
  model_factor        *= 1.0 / (4*M_PI*pow(AU_2_GEV, 2));
  
  // converting to ( 1 / m^2 s)
  model_factor        *= pow(M_2_GEV, 2)*SEC_2_GEV ;

  return model_factor ;
}

//--------------------------------------------------------------
// Function for plotting contours (inherited from base class)
// Threshold case in (m_A', m_DM) plane
double DarkPhoton::ContFuncThresh_LM(double m_dp, double m_dm)
{
  PROFILE_FUNCTION() ;


  if(GetDecayProd().size()*active_bin.GetERange().min > m_dm ||
     GetDecayProd().size()*active_bin.GetERange().max < m_dm) return 1000 ;

  dark_model.eps = 1e-9 ;

  dark_model.mDM = m_dm ;
  dark_model.mDP = m_dp ;

  // Just for the flag mostly
  // SetModelPars(dark_model) ;

  FixRelic();

  double  model_factor = 1.0 ; // Br = 1
  model_factor        *= GetDecayProbT()     ;
  // model_factor        *= GetDecayProdNum()   ;
  
  // ( X X ==> fi fi ) Annihilation
  model_factor        *= 2*GetAnnRate()        ;
  model_factor        *= 1.0 / (4*M_PI*pow(AU_2_GEV, 2));
  
  // converting to ( 1 / m^2 s)
  model_factor        *= pow(M_2_GEV, 2)*SEC_2_GEV ;

  return model_factor ;
}

//--------------------------------------------------------------
// Function for plotting contours (inherited from base class)
// Threshold case in (m_DM, eps) plane
double DarkPhoton::ContFuncThresh_MG(double m_dm, double in_eps)
{
  PROFILE_FUNCTION() ;

  if(GetDecayProd().size()*active_bin.GetERange().min > m_dm ||
     GetDecayProd().size()*active_bin.GetERange().max < m_dm) return 1000 ;

  dark_model.eps = in_eps ;
  dark_model.mDM = m_dm ;

  dark_model.mDP = 4*ELECTRON_M_GEV ;

  // Just for the flag mostly
  // SetModelPars(dark_model) ;

  FixRelic();

  double  model_factor = 1.0 ; // Br = 1
  model_factor        *= GetDecayProbT()     ;
  // model_factor        *= GetDecayProdNum()   ;
  
  // ( X X ==> fi fi ) Annihilation
  model_factor        *= 2*GetAnnRate()        ;
  model_factor        *= 1.0 / (4*M_PI*pow(AU_2_GEV, 2));
  
  // converting to ( 1 / m^2 s)
  model_factor        *= pow(M_2_GEV, 2)*SEC_2_GEV ;

  return model_factor ;
}

//--------------------------------------------------------------
// Function for plotting contours (inherited from base class)
// Boosted case in (m_A', eps) plane
double DarkPhoton::ContFuncBoost_LG(double x, double y)
{
  PROFILE_FUNCTION() ;
  dark_model.eps = y ;
  dark_model.mDP = x ;

  // Just for the flag mostly
  // SetModelPars(dark_model) ;

  FixRelic();

  double  model_factor = GetDPhotonBr()       ;
  model_factor        *= GetDecayProbB()      ;
  // model_factor        *= GetDecayProdNum()    ;
  model_factor        *= 2*GetAnnRate()         ;
  model_factor        *= 1.0 / (4*M_PI*pow(AU_2_GEV, 2)) ;
  
  // converting to ( 1 / m^2 s)
  model_factor        *= pow(M_2_GEV, 2)*SEC_2_GEV ; 

  return model_factor ;

}


//--------------------------------------------------------------
// Function for plotting contours (inherited from base class)
// Boosted case in (m_A', m_dm) plane
// !Not implemented yet!
double DarkPhoton::ContFuncBoost_LM(double x, double y)
{
  // PROFILE_FUNCTION() ;
  // dark_model.eps = 1e-9 ;

  // dark_model.mDM = y ;
  // dark_model.mDP = x ;

  // // Just for the flag mostly
  // // SetModelPars(dark_model) ;

  // FixRelic();

  // double  model_factor = GetDPhotonBr()       ;
  // model_factor        *= GetDecayProbB()      ;
  // // model_factor        *= GetDecayProdNum()    ;
  // model_factor        *= 2*GetAnnRate()         ;
  // model_factor        *= 1.0 / (4*M_PI*pow(AU_2_GEV, 2)) ;
  
  // // converting to ( 1 / m^2 s)
  // model_factor        *= pow(M_2_GEV, 2)*SEC_2_GEV ; 

  // return model_factor ;
  return 1 ;
}

//--------------------------------------------------------------
// Function for plotting contours (inherited from base class)
// Boosted case in (m_dm, eps) plane
// !Not implemented yet!
double DarkPhoton::ContFuncBoost_MG(double x, double y)
{
  // PROFILE_FUNCTION() ;

  // dark_model.eps = y ;
  // dark_model.mDM = x ;

  // dark_model.mDP = 4*ELECTRON_M_GEV ;

  // // Just for the flag mostly
  // // SetModelPars(dark_model) ;

  // FixRelic();

  // double  model_factor = GetDPhotonBr()       ;
  // model_factor        *= GetDecayProbB()      ;
  // // model_factor        *= GetDecayProdNum()    ;
  // model_factor        *= 2*GetAnnRate()         ;
  // model_factor        *= 1.0 / (4*M_PI*pow(AU_2_GEV, 2)) ;
  
  // // converting to ( 1 / m^2 s)
  // model_factor        *= pow(M_2_GEV, 2)*SEC_2_GEV ; 

  // return model_factor ;
  return 1;
}

//--------------------------------------------------------------
// Initializing the model
void DarkPhoton::Init()
{
  // DPhotPar dphot_params ;
  // Excact value doesn't matter
  // dphot_params.eps = 1e-10 ;
  // dphot_params.mDM = 2000  ;
  // dphot_params.mDP = 0.1   ;

  SetDecayProd({DMSolarSignal::Electron, DMSolarSignal::Positron}) ;
  // SetModelPars(dphot_params) ;
  // FixRelic();

  Zaki::String::Directory dir(__FILE__) ;
  if(dir.ParentDir().ParentDir().Str() == dir.ParentDir().Str())
    dir = "" ;
  else
    dir = dir.ParentDir().ParentDir() + "/" ;

  ImportDarkPhotonBr(dir + "data/DPhoton_Br.dat");
}

//--------------------------------------------------------------
// Overriding the method from the base class
void DarkPhoton::SetDMMass(const double in_mass) 
{
  dark_model.mDM = in_mass;
  Model::SetDMMass(in_mass) ;

  // band-aid! needs to be changed! 
  if(DM_mass_set.size() == 0)
    DM_mass_set.push_back(in_mass) ;
}

//--------------------------------------------------------------
// Overriding the method from the base class
void DarkPhoton::SetDMMass(const std::vector<double>& in_mass) 
{
  Model::SetDMMass(in_mass) ;
}
//==============================================================
