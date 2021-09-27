/*
  Model class

*/

// ROOT
// #include <TF1.h>
// #include <TCanvas.h>
// #include <TStyle.h>
// #include <Math/WrappedFunction.h>
// #include <Math/Integrator.h>
// #include <Math/Interpolator.h>

// Local headers
#include "DMSS/Model.hpp"

//==============================================================
//--------------------------------------------------------------
/// Destructor
Model::~Model() { }

//--------------------------------------------------------------
/// Default constructor
Model::Model()  
  : Prog("Model", true)  
{
  // logger.SetUnit("Model") ;
  // SetName("[NOT SET]") ;

  /// Default plot options
  // plt_opt = {1000, 1000, "X", "Y", false} ;
}

//--------------------------------------------------------------
/// Constructor via name
Model::Model(const std::string& in_name)
{ 
  SetName(in_name) ;
  
  /// Default plot options
  // plt_opt = {1000, 1000, "X", "Y", false} ;
}

//--------------------------------------------------------------
// Prints the model info
// void Model::Print() const
// {
//   std::cout << "* ---------------------------------------- *" << "\n";
//   std::cout << "            Model:   "<< GetName() << "\n";
//   // for (size_t i = 0; i < params.size() ; i++)
//   // {
//   //   std::cout << params[i] << "\n";
//   // }
//   std::cout << "* ---------------------------------------- *" << "\n";
// }

//--------------------------------------------------------------
/// Init in the base class
void Model::Init() {}

//--------------------------------------------------------------
void Model::SetTGrid(const PlotMode& in_mode, 
                     const Zaki::Math::Grid2D& in_grid, 
                     const PlotOptions& in_pltOpts,
                     const bool in_fixed) 
{
  thresh_grids[in_mode] = {in_grid, in_pltOpts, in_fixed} ;
  // set_thresh_grid_flag = true ;
}

//--------------------------------------------------------------
void Model::SetBGrid(const PlotMode& in_mode,
                     const Zaki::Math::Grid2D& in_grid,
                     const PlotOptions& in_pltOpts,
                     const bool in_fixed) 
{
  boos_grids[in_mode] = {in_grid, in_pltOpts, in_fixed} ;
  // set_boos_grid_flag = true ;
}

//--------------------------------------------------------------
PlotGrid* Model::GetBGrid(const PlotMode& in_mode)
{
  if(boos_grids.count(in_mode))
  {
    return &boos_grids[in_mode];
  }
  else
  {
    Z_LOG_ERROR("Boosted grid is not set, use 'Model::SetBGrid(const PlotMode&)' first!") ;
    return nullptr ;
  }

}

//--------------------------------------------------------------
PlotGrid* Model::GetTGrid(const PlotMode& in_mode)
{
  if (thresh_grids.count(in_mode))
  {
    return &thresh_grids[in_mode];
  }
  else
  {
    Z_LOG_ERROR("Threshold grid is not set, use 'Model::SetTGrid(const PlotMode&)' first!") ;
    return nullptr ;
  }
}

//--------------------------------------------------------------
void Model::SetDecayProd(const std::vector<DMSolarSignal::Particle>& in)
{
  decay_products = in ;
}

//--------------------------------------------------------------
// Number of decay products of the mediator (fi)
// DM + DM => fi + fi ; fi => e+e- ...
std::vector<DMSolarSignal::Particle> Model::GetDecayProd() const
{
  return decay_products ;
}

//--------------------------------------------------------------
void Model::SetECut(double e_cut_in)
{
  energy_cut = e_cut_in ;
}

//--------------------------------------------------------------
double Model::GetECut() const
{
  return energy_cut;
}

//--------------------------------------------------------------
// Setting the dark matter mass
void Model::SetDMMass(const double in_mass) 
{
  mDM = in_mass;

  // band-aid fixing, needs to be changed! 
  if(DM_mass_set.size() == 0)
    DM_mass_set.push_back(in_mass) ;
}

//--------------------------------------------------------------
// Setting the DM mass for multiple contours in the boosted case
void Model::SetDMMass(const std::vector<double>& in_mset) 
{
  mDM = in_mset[0] ;
  DM_mass_set = in_mset ;
}

//--------------------------------------------------------------
// Setting the dark matter mass
double Model::GetDMMass() const
{
  return mDM ;
}

//--------------------------------------------------------------
// Getting the dark matter mass Set
std::vector<double> Model::GetDMMassSet() const
{
  return DM_mass_set ;
}
    

//--------------------------------------------------------------
// The spectrum of the model in a given bin
double Model::GetSpectrum(Bin& in_bin) 
{
  // Returns the starting energy of the bin
  // if greater than the energy cut
  // else it returns the energy cut
  double e1 = std::max({in_bin.GetERange().min, energy_cut}) ;
  double e2 = std::min({in_bin.GetERange().max, mDM}) ;

  return (e2 - e1) / mDM ;
}

//--------------------------------------------------------------
/// Sets the bin that is under consideration (active)
void Model::SetActiveBin(const Bin& in_bin) 
{
  active_bin = in_bin ;
}

//--------------------------------------------------------------
/// Gets the bin that is under consideration (active)
Bin Model::GetActiveBin() const 
{
  return active_bin ;
}

//--------------------------------------------------------------
// Setting plot options
// void Model::SetPlotOptions(const PlotOptions& in_plot_opt) 
// {
//   plt_opt = in_plot_opt ;
//   set_plt_opts_flag = true ;

// }

//--------------------------------------------------------------
// Returns plot options
// PlotOptions Model::GetPlotOptions() 
// {
//   if (!set_plt_opts_flag)
//     Z_LOG_WARNING("Plot options are not set for this model, using the default options instead.") ;

//   return plt_opt ;
// }

//--------------------------------------------------------------
// Checking if the grid values are fixed with respect to changes
// in another parameter (e.g. M_DM for (L_dec, G_ann) plane)
// bool Model::HasFixedGrid(const PlotMode& in_mode) 
// {
//   return has_fixed_grid ;
// }

//--------------------------------------------------------------
// Fixing the grid so that it returns true if the grid values 
// are fixed with respect to the DM mass
// void Model::SetFixedGrid(bool in_bool) 
// {
//   has_fixed_grid = in_bool ;
// }

//==============================================================
