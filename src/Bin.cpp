/*
  Bin class

*/

#include <TCanvas.h>
#include <TStyle.h>

#include "DMSS/Bin.hpp"

//==============================================================
std::ostream& operator << ( std::ostream &output, Flux f )
{ 
  char tmp[150] ;
  sprintf(tmp, "%.3e \u00B1 %.3e(stat) \u00B1 %.3e(sys)", f.val, f.statErr, f.sysErr) ;
  output << tmp ;
  return output;            
}

//==============================================================

//--------------------------------------------------------------
/// Constructor 1
Bin::Bin()  
  : Prog("Bin", true)
{ }

//--------------------------------------------------------------
/// Full Constructor
Bin::Bin(const Zaki::Math::Range<double>& in_e_range,
    const Zaki::Math::Quantity& in_e_center,
    const Zaki::Math::Quantity& in_n_count,
    const Flux& in_flux)
    : Prog("Bin", true),
    flux(in_flux), e_range(in_e_range),
    e_center(in_e_center), n_count(in_n_count),
    set_flux_flag(true), set_e_range_flag(true),
    set_e_center_flag(true), set_n_count_flag(true)
{} 

//--------------------------------------------------------------
/// Copy constructor
Bin::Bin(const Bin& other)
    : Prog("Bin", true),
    flux(other.flux),
    e_range(other.e_range),
    e_center(other.e_center),
    n_count(other.n_count),
    set_flux_flag(other.set_flux_flag),
    set_e_range_flag(other.set_e_range_flag),
    set_e_center_flag(other.set_e_center_flag),
    set_n_count_flag(other.set_n_count_flag),
    divided_bin_flag(other.divided_bin_flag),
    t_bin_obs_set(other.t_bin_obs_set),
    t_bin_chops(other.t_bin_chops)
{
    t_bin_hist = *(static_cast<TH1F*>(other.t_bin_hist.Clone())) ;
    t_bin_hist.SetName(("t_bin_hist_" + std::to_string(rand() % 990000)).c_str()) ;
//    Z_LOG_INFO("Bin copy constructor called -> from " + std::string(other.t_bin_hist.GetName()) + " to " + std::string(t_bin_hist.GetName()) ) ;
}

//--------------------------------------------------------------
Bin::~Bin() {} 

//--------------------------------------------------------------
// [] Operator
// std::string const& CSVRow::operator[](std::size_t index) const 
// {
//   return m_data[index];
// }


//--------------------------------------------------------------
// size method

//--------------------------------------------------------------
// void Bin::SetLabel(const std::string& label_in) 
// {
//   set_label_flag = true ;
//   label = label_in ;
// }

//--------------------------------------------------------------
// SetFlux method 
void Bin::SetFlux(double fi, double stat_err, double sys_err) 
{
  set_flux_flag = true ;
  flux = {fi, stat_err, sys_err} ; 
}
//--------------------------------------------------------------
void Bin::SetERange(double E_min, double E_max) 
{
  set_e_range_flag = true ;
  // e_range = std::make_pair(E_min, E_max); 
  e_range = {E_min, E_max} ;
}
//--------------------------------------------------------------
void Bin::SetECenter(double Ec, double Ec_err) 
{
  set_e_center_flag = true ;
  // e_center = std::make_pair(Ec, Ec_err); 
  e_center = {Ec, Ec_err} ;
}
//--------------------------------------------------------------
void Bin::SetNCount(double n, double n_err) 
{
  set_n_count_flag = true ;
  // n_count = std::make_pair(n, n_err); 
  n_count = {n , n_err} ;
}
//--------------------------------------------------------------
// std::string Bin::GetLabel() 
// {
//   if (set_label_flag)
//     return label ;
//   else 
//     return "[NOT SET]" ;
// }
//--------------------------------------------------------------
Flux Bin::GetFlux() const
{
  return flux ;
}
//--------------------------------------------------------------
Zaki::Math::Range<double> Bin::GetERange() const
{
  return e_range ;
}
//--------------------------------------------------------------
Zaki::Math::Quantity Bin::GetECenter() const
{
  return e_center ;
}
//--------------------------------------------------------------
Zaki::Math::Quantity Bin::GetNCount() const
{
  return n_count ;
}
//--------------------------------------------------------------
void Bin::Print() const
{
  std::cout<<"\n ------------------------------------------------------------------------\n" ;
  std::cout<<"| Bin Label: "<<GetName() << "\n";
  std::cout<<"| Energy Range: "<< e_range << " GeV"<< "\n";
  std::cout<<"| Characteristic Energy: " << e_center <<" GeV"<< "\n";
  std::cout<<"| Counts: "<< n_count << "\n";
  std::cout<<"| Flux: "<< flux << "  [m^2 sr s GeV]^−1"<< "\n";
  std::cout<<" ------------------------------------------------------------------------\n" ;
}

//--------------------------------------------------------------
/// String form of the bin (for printing as a table)
std::string Bin::Str() const 
{
  char tmp[300] ;
  sprintf(tmp, "%6.1f, %6.1f, %6.1f, %4.1f, %5.0f, %3.0f, %5.3e, %5.3e, %5.3e",
   e_range.min, e_range.max, e_center.val, e_center.err, 
   n_count.val, n_count.err, flux.val, flux.statErr, flux.sysErr ) ;

  return std::string(tmp) ;      
}

//--------------------------------------------------------------
size_t Bin::GetTBinChops() const
{
  return t_bin_chops;
}
//--------------------------------------------------------------
TH1F Bin::divide(size_t chops)
{
  if (!set_n_count_flag) 
  { 
    Z_LOG_ERROR("Number of counts has not been set yet!") ;
  }

  t_bin_chops = chops ;

  // size_t bin_num        = chops ;
  double frac             = 1.0 / chops ;
  size_t t_min            = 0           ; 
  size_t t_max            = chops       ; 

  Zaki::Math::Quantity tmp_obs = {frac*n_count.val, pow(frac, 0.5)*n_count.err } ;

  // std::cout<<"\n -Bin.Divide():-->>> Obs. Val: "<<obs_val<<", Obs. Err: "<<obs_err<<", Bin #: "<<chops<<" \n";

  char title_char[200] ;    
  sprintf(title_char, "Constant Observed Events (E = %.1f GeV)", e_center.val) ;

    std::string tmp_str = "obs_" + std::to_string(rand() % 990000) ;
    TH1F o(tmp_str.c_str(), title_char, chops, t_min,t_max);
//    Z_LOG_NOTE("Created TH1F: " + tmp_str) ;
    
  t_bin_obs_set     = {} ;
  t_bin_obs_set.reserve(chops) ;

  // Filling the observed histogram
  for(size_t i=0; i < chops; i++)
	{
    o.SetBinContent(i+1, tmp_obs.val) ;
    o.SetBinError(i+1, tmp_obs.err) ;
    t_bin_obs_set.push_back(tmp_obs);
  }

  t_bin_hist = o ;

  divided_bin_flag = true ; 

  return o;
}
//--------------------------------------------------------------
TH1F Bin::GetTBinHist() const
{
  if (!divided_bin_flag )
    Z_LOG_ERROR("Bin hasn't been divided yet, use 'divide(double chops)' first!") ;
  
  return t_bin_hist;
}
//--------------------------------------------------------------
std::vector<Zaki::Math::Quantity> Bin::GetTBinObsSet() const
{
  return t_bin_obs_set;
}

//--------------------------------------------------------------
void Bin::Plot() const
{
  // TH1F h = this->divide(chops);
  TH1F h = GetTBinHist() ;

  h.GetXaxis()->SetTitle("Cycle");
  h.GetYaxis()->SetTitle("Observed Counts");

  h.SetLineWidth(2) ;

  TCanvas c("c", GetName().c_str(), 1000, 600) ;
  c.SetGrid();
  
  // c1->Divide(2) ;
  // c1->cd(1) ;

  gPad->SetGrid() ;
  gStyle->SetOptStat(0);

  h.Draw("E1") ;

  // Output file
  char out_file_char[200] ;    
  sprintf(out_file_char, "E_%.1f_GeV_%d_Cycles.pdf", e_center.val, static_cast<int>(t_bin_chops)) ;

  c.SaveAs((wrk_dir + "/" + out_file_char).Str().c_str()) ;
}
//==============================================================
