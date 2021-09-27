/*
  Data class

*/

#include <stdlib.h>

// Root
#include <TH1F.h>
#include <TCanvas.h>
#include <TStyle.h>

#include <Zaki/File/CSVIterator.hpp>
#include <Zaki/File/VecSaver.hpp>

// Local headers
#include "DMSS/Data.hpp"

//==============================================================

//--------------------------------------------------------------
// Constructor
Data::Data()   : Prog("Data", true) {}

//--------------------------------------------------------------
// [] Operator
Bin& Data::operator[](size_t index)
{
  return bins[index];
}

//--------------------------------------------------------------
// size
size_t Data::size() const 
{
  return bins.size() ;
}

//--------------------------------------------------------------
// Adding bins
void Data::AddBin(const Bin& b) 
{
  bins.push_back(b) ;

  // For the first bin we add the left edge
  if(bins.size() == 1)
    edge_vec.push_back(b.GetERange().min);

  // For all bins we add the right edge
  edge_vec.push_back(b.GetERange().max);
}

//--------------------------------------------------------------
// Adding bins
void Data::AddBin(const Zaki::Math::Range<double>& in_e_range,
            const Zaki::Math::Quantity& in_e_center,
            const Zaki::Math::Quantity& in_n_count,
            const Flux& in_flux)
{
    bins.emplace_back(in_e_range, in_e_center, in_n_count, in_flux) ;
    
    // For the first bin we add the left edge
    if(bins.size() == 1)
        edge_vec.emplace_back(in_e_range.min);
    
    // For all bins we add the right edge
    edge_vec.emplace_back(in_e_range.max);
}

//--------------------------------------------------------------
/// Scales the number of events (and flux) by r,
/// Scales errors by sqrt(r)
void Data::ScaleData(const double& r ) 
{
  for (auto &&b : bins)
  {
    b.SetNCount(b.GetNCount().val*r, b.GetNCount().err*sqrt(r)) ;
    b.SetFlux(b.GetFlux().val*r, 
      b.GetFlux().statErr*sqrt(r), b.GetFlux().sysErr*sqrt(r)) ;
  }
  
}

//--------------------------------------------------------------
void Data::Input(const Zaki::String::Directory& file_name)
{  
  std::ifstream     file(file_name.Str());
  
  if (file.fail())
  {
    Z_LOG_ERROR(("Data file '" + file_name.Str() + "' cannot be opened!").c_str()) ;
    exit(EXIT_FAILURE) ;
    // return;
  }
  else
  {
    Z_LOG_INFO(("Data file '" + file_name.Str() + "' imported.").c_str()) ;
  }
  
  // Reading the input file
  for(Zaki::File::CSVIterator loop(file); loop != Zaki::File::CSVIterator(); ++loop)
  {
//    Bin b0 ;
//    // b0.SetName(GetName().c_str()) ;
//    b0.SetERange(std::stof((*loop)[0]), std::stof((*loop)[1])) ;
//    b0.SetECenter(std::stof((*loop)[2]), std::stof((*loop)[3])) ;
//    b0.SetNCount(std::stof((*loop)[4]), std::stof((*loop)[5])) ;
//    b0.SetFlux(std::stof((*loop)[6]), std::stof((*loop)[7]), std::stof((*loop)[8]) ) ;
      
    AddBin(
            { std::stof((*loop)[0]), std::stof((*loop)[1]) },
            { std::stof((*loop)[2]), std::stof((*loop)[3]) },
            { std::stof((*loop)[4]), std::stof((*loop)[5]) },
            { std::stof((*loop)[6]), std::stof((*loop)[7]), std::stof((*loop)[8]) }
           ) ;
//    AddBin(b0);
  }
}

//--------------------------------------------------------------
void Data::SetStates(const std::vector<DMSolarSignal::Particle>& in_states)
{
  states = in_states;
}

//--------------------------------------------------------------
std::vector<DMSolarSignal::Particle> Data::GetStates() const
{
  return states ;
}

//--------------------------------------------------------------
void Data::Print() const
{
  std::cout<<"\n ------------------------------------------------\n" ;
  std::cout<<"| Data Label: "<<GetName() << "\n";
  std::cout<<"| Total Bins: "<<bins.size()<< "\n";
  std::cout<<" ------------------------------------------------\n" ;

  // if(opt =="a")
  PrintBins(bins.size()) ;
}

//--------------------------------------------------------------
void Data::PrintBins(size_t n) const
{
  if ( n > bins.size() )
  {
    Z_LOG_WARNING("Number of bins requested larger than the data size, printing '"
                  + std::to_string(bins.size()) + "' bins instead...") ;
  }
  for(size_t i=0 ; i<bins.size() && i<n ; ++i)
  {
    std::cout<<"\n Bin # "<<i+1<<"";
    bins[i].Print() ;
  }
}

//--------------------------------------------------------------
void Data::Plot(const Zaki::String::Directory& f_name) const
{
  double edge_arr[bins.size() + 1] ;

  std::copy(edge_vec.begin(), edge_vec.end(), edge_arr);

  TH1F o("observed", "Total Counts", bins.size(), edge_arr);

  // Filling the observed histogram
  for(size_t i=0; i<bins.size(); i++)
	{
    o.SetBinContent(i+1, bins[i].GetNCount().val) ;
    o.SetBinError(i+1,  bins[i].GetNCount().err) ;
  }

  o.GetXaxis()->SetTitle("E  [GeV]");
  o.GetYaxis()->SetTitle("Observed Counts");

  o.SetLineWidth(2) ;

  TCanvas c("c", GetName().c_str(), 800, 600) ;
  c.SetGrid();
  
  // c1->Divide(2) ;
  // c1->cd(1) ;

  c.SetLogy();

  // gPad->SetGrid() ;
  gStyle->SetOptStat(0);

  o.Draw("E1") ;

  // Output file
  c.SaveAs((wrk_dir + f_name + ".pdf").Str().c_str()) ;
    
}

//--------------------------------------------------------------
std::vector<Bin> Data::GetBins() const
{
  return bins;
}

//--------------------------------------------------------------
void Data::Export(const Zaki::String::Directory& dir) 
{
  Zaki::File::VecSaver my_saver(wrk_dir + dir) ;
  my_saver.Export1D(bins) ;

  Z_LOG_INFO("Data exported to '"+(wrk_dir + dir).Str()+"'.") ;
}
//--------------------------------------------------------------

//==============================================================
