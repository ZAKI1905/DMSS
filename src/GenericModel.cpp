/*
  GenericModel class

*/

#include <Zaki/Physics/Constants.hpp>

// Local headers
#include "DMSS/GenericModel.hpp"

//==============================================================
using namespace Zaki::Physics ;
//--------------------------------------------------------------
// Constructor
GenericModel::GenericModel() 
{
  size_t tmp_Res = 200; 
  // logger.SetUnit("GenericModel") ;
  SetName("Generic") ;

  //....................................................................//
  //                  Setting the boosted grids                         //
  //....................................................................//
  SetBGrid(PlotMode::Ldec_Gann , 
           {{{5e4, 1e12}, tmp_Res, "Log"}, {{1e-5, 9e-1}, tmp_Res, "Log"}},
           {1000, 1000, "L_{dec} [km]", "#Gamma [GeV]", true},
           false) ;
  SetBGrid(PlotMode::Ldec_Mdm , 
           {{{5e4, 1e12}, tmp_Res, "Log"}, {{100, 4000}, tmp_Res, "Linear"}},
           {1000, 1000, "L_{dec} [km]", "M_{DM} [GeV]", true},
           false) ;
  SetBGrid(PlotMode::Mdm_Gann , 
            {{{100, 4000}, tmp_Res, "Linear"}, {{3e-5, 9e-1}, tmp_Res, "Log"}},
            {1000, 1000, "M_{DM} [GeV]", "#Gamma [GeV]", true},
            false) ;
  //....................................................................//

  //....................................................................//
  //                  Setting the threshold grids                       //
  //....................................................................//
  SetTGrid(PlotMode::Ldec_Gann,
           {{{5e4, 1e12}, tmp_Res, "Log"}, {{1e-6, 9e-1}, tmp_Res, "Log"}},
           {1000, 1000, "L_{dec} [km]", "#Gamma [GeV]", true},
          true) ;
  SetTGrid(PlotMode::Ldec_Mdm , 
           {{{5e4, 1e12}, tmp_Res, "Log"}, {{100, 4000}, tmp_Res, "Log"}},
           {1000, 1000, "L_{dec} [km]", "M_{DM} [GeV]", true},
           false) ;
  SetTGrid(PlotMode::Mdm_Gann , 
            {{{100, 4000}, tmp_Res, "Log"}, {{3e-6, 9e-1}, tmp_Res, "Log"}},
            {1000, 1000, "M_{DM} [GeV]", "#Gamma [GeV]", true},
            false) ;
  //....................................................................//

  // Setting the plot options:
  //  {width, height, x-Label, y-label, connected}
  // SetPlotOptions({1000, 1000, "L_{dec} [km]", "#Gamma [GeV]", false}) ;

  // Setting the decay products of the mediator
  SetDecayProd({DMSolarSignal::Electron, DMSolarSignal::Positron}) ;

  // // Fixing the grid to save time in threshold analysis
  // SetFixedGrid() ;

  Z_LOG_NOTE(("GenericModel constructor called for: " + PtrStr()).c_str()) ;

}

//--------------------------------------------------------------
// IClone method
GenericModel* GenericModel::IClone() const
{

  Z_LOG_NOTE(("GenericModel::Clone() called for: " + PtrStr()).c_str()) ;

  return (new GenericModel(*this));
}

//--------------------------------------------------------------
// Destructor
GenericModel::~GenericModel() 
{
  Z_LOG_NOTE(("GenericModel destructor called for: " + PtrStr()).c_str()) ;
}

//--------------------------------------------------------------
// Function for plotting contours (inherited from base class)
// Threshold case in (L_dec, G_ann) plane
double GenericModel::ContFuncThresh_LG(double dec_len, double gamma)
{
  PROFILE_FUNCTION() ;

  // Changing units from km to GeV
  double L_dec     = dec_len * KM_2_GEV ;

  // Decay probability in the threshold case
  double p  = exp( - SUN_R_GEV / L_dec ) ;
  p        -= exp( - EARTH_2_SUN_GEV / L_dec ) ;

  double  model_factor = 2*gamma*p ;
  model_factor        *= 1.0 / (4*M_PI*pow(AU_2_GEV, 2)) ;
  
  // converting to ( 1 / m^2 s)
  model_factor        *= pow(M_2_GEV, 2)*SEC_2_GEV ;

  return model_factor ;
}

//--------------------------------------------------------------
// Function for plotting contours (inherited from base class)
// Threshold case in (M_DM, G_ann) plane
double GenericModel::ContFuncThresh_MG(double m_dm, double gamma)
{
  PROFILE_FUNCTION() ;

  if(GetDecayProd().size()*active_bin.GetERange().min > m_dm ||
     GetDecayProd().size()*active_bin.GetERange().max < m_dm) return 1000 ;

  double dec_len = 2e07;

  // Changing units from km to GeV
  double L_dec     = dec_len * KM_2_GEV ;

  // Decay probability in the threshold case
  double p  = exp( - SUN_R_GEV / L_dec ) ;
  p        -= exp( - EARTH_2_SUN_GEV / L_dec ) ;

  double  model_factor = 2*gamma*p ;
  model_factor        *= 1.0 / (4*M_PI*pow(AU_2_GEV, 2)) ;
  
  // converting to ( 1 / m^2 s)
  model_factor        *= pow(M_2_GEV, 2)*SEC_2_GEV ;

  return model_factor ;
}

//--------------------------------------------------------------
// Function for plotting contours (inherited from base class)
// Threshold case in (L_dec, M_DM) plane
double GenericModel::ContFuncThresh_LM(double dec_len, double m_dm)
{
  PROFILE_FUNCTION() ;

  if(GetDecayProd().size()*active_bin.GetERange().min > m_dm ||
     GetDecayProd().size()*active_bin.GetERange().max < m_dm) return 1000 ;

  double gamma = 1e-3;

  // Changing units from km to GeV
  double L_dec     = dec_len * KM_2_GEV ;

  // Decay probability in the threshold case
  double p  = exp( - SUN_R_GEV / L_dec ) ;
  p        -= exp( - EARTH_2_SUN_GEV / L_dec ) ;

  double  model_factor = 2*gamma*p ;
  model_factor        *= 1.0 / (4*M_PI*pow(AU_2_GEV, 2)) ;
  
  // converting to ( 1 / m^2 s)
  model_factor        *= pow(M_2_GEV, 2)*SEC_2_GEV ;

  return model_factor ;
}

//--------------------------------------------------------------
// Function for plotting contours (inherited from base class)
// Boosted case in (L_dec, G_ann) plane
double GenericModel::ContFuncBoost_LG(double dec_len, double gamma)
{
  PROFILE_FUNCTION() ;
  // Changing units from km to GeV
  double L_dec     = dec_len * KM_2_GEV ;

  // Decay probability in the threshold case
  double p  = exp( - SUN_R_GEV / L_dec )       ;
  p        -= exp( - EARTH_2_SUN_GEV / L_dec ) ;

  double  model_factor = 2*gamma*p      ;
  model_factor        *= 1.0 / (4*M_PI*pow(AU_2_GEV, 2)) ;
  
  // converting to ( 1 / m^2 s)
  model_factor        *= pow(M_2_GEV, 2)*SEC_2_GEV ; 

  return model_factor ;
}

//--------------------------------------------------------------
// Function for plotting contours (inherited from base class)
// Boosted case in (M_DM, G_ann) plane
// ! Not implemented yet !
double GenericModel::ContFuncBoost_MG(double m_dm, double gamma)
{
  // PROFILE_FUNCTION() ;

  // double dec_len = 2e07;

  // // Changing units from km to GeV
  // double L_dec     = dec_len * KM_2_GEV ;

  // // Decay probability in the threshold case
  // double p  = exp( - SUN_R_GEV / L_dec )       ;
  // p        -= exp( - EARTH_2_SUN_GEV / L_dec ) ;

  // double  model_factor = 2*gamma*p      ;
  // model_factor        *= 1.0 / (4*M_PI*pow(AU_2_GEV, 2)) ;
  
  // // converting to ( 1 / m^2 s)
  // model_factor        *= pow(M_2_GEV, 2)*SEC_2_GEV ; 

  // return model_factor ;

  return 1;
}

//--------------------------------------------------------------
// Function for plotting contours (inherited from base class)
// Boosted case in (L_dec, M_dm) plane
// ! Not implemented yet !
double GenericModel::ContFuncBoost_LM(double dec_len, double m_dm)
{
  // PROFILE_FUNCTION() ;

  // double gamma = 1e-3;

  // // Changing units from km to GeV
  // double L_dec     = dec_len * KM_2_GEV ;

  // // Decay probability in the threshold case
  // double p  = exp( - SUN_R_GEV / L_dec )       ;
  // p        -= exp( - EARTH_2_SUN_GEV / L_dec ) ;

  // double  model_factor = 2*gamma*p      ;
  // model_factor        *= 1.0 / (4*M_PI*pow(AU_2_GEV, 2)) ;
  
  // // converting to ( 1 / m^2 s)
  // model_factor        *= pow(M_2_GEV, 2)*SEC_2_GEV ; 

  // return model_factor ;
  return 1;
}

//==============================================================
