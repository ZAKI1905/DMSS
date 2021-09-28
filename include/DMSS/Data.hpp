#ifndef DMSS_Data_H
#define DMSS_Data_H

#include "DMSS/Bin.hpp"

//==============================================================
namespace DMSolarSignal{
enum Particle
{
  Electron = 0, Positron, Proton 
} ;
}
//==============================================================
class Data : public Prog
{
  //--------------------------------------------------------------
  public:
    Data() ;

  //............................................
  // Setters
  //............................................
    // Adding an individual bin
    void AddBin(const Bin&) ;
    void AddBin(const Zaki::Math::Range<double>& in_e_range,
                const Zaki::Math::Quantity& in_e_center,
                const Zaki::Math::Quantity& in_n_count,
                const Flux& in_flux);

    // Adding a set of bins from file
    void Input(const Zaki::String::Directory&) ;

    // Setting the states
    void SetStates(const std::vector<DMSolarSignal::Particle>&) ;
  //............................................

  //............................................
  // Getters
  //............................................
    std::vector<Bin> GetBins() const ;
    std::vector<DMSolarSignal::Particle> GetStates() const ;

    // Printing the data set info
    void Print() const override;

    // Printing the first 'n' bins info
    void PrintBins(size_t n = 3) const ;

    void Plot(const Zaki::String::Directory&) const ;
    Bin& operator[](size_t) ;
    size_t size() const ;

    void ScaleData(const double& ) ;
    void Export(const Zaki::String::Directory& dir) ;

  //--------------------------------------------------------------
  private:
    std::vector<Bin> bins ;
    std::vector<double> edge_vec ;
    std::vector<DMSolarSignal::Particle> states;
};

//==============================================================
#endif /*DMSS_Data_H*/
