#ifndef DMSS_Bin_H
#define DMSS_Bin_H

// Root
#include <TH1F.h>

// Zaki
#include <Zaki/Math/Math_Core.hpp>

// Local headers
#include "DMSS/Prog.hpp"

//==============================================================
struct Flux
{
  double val, statErr, sysErr ; 
};

std::ostream& operator << ( std::ostream &output, Flux f );
//==============================================================
class Bin : public Prog
{

  //--------------------------------------------------------------
  public:

    // Constructor 1
    Bin() ;
    
    // Constructor 2
    Bin(const Zaki::Math::Range<double>& in_e_range,
        const Zaki::Math::Quantity& in_e_center,
        const Zaki::Math::Quantity& in_n_count,
        const Flux& in_flux) ;
    
    Bin(const Bin& other) ;
    
    ~Bin() ;

    //............................................
    // Setters
    //............................................
    // void SetLabel(const std::string&) ;
    void SetFlux(double,double,double) ;
    void SetERange(double, double) ;
    void SetECenter(double, double) ;
    void SetNCount(double, double) ;
    //............................................

    //............................................
    // Getting the variables
    //............................................
    Flux GetFlux() const;
    Zaki::Math::Range<double> GetERange() const;
    Zaki::Math::Quantity GetECenter() const;
    Zaki::Math::Quantity GetNCount() const;
    std::vector<Zaki::Math::Quantity> GetTBinObsSet() const;
    // std::vector<double> GetTBinObsSetErr() const;
    TH1F GetTBinHist() const;
    size_t GetTBinChops() const;
    //............................................

    // Dividing one bin into more bins
    TH1F divide(size_t) ;

    //............................................
    // Plotters
    //............................................
    // Printing the bin info
    void Print() const override ;

    /// String form of the bin
    std::string Str() const ;

    // Plot the time-binned result
    void Plot() const ;

  //--------------------------------------------------------------
  private:

    Flux flux;
    Zaki::Math::Range<double> e_range  ;
    Zaki::Math::Quantity e_center ;
    Zaki::Math::Quantity n_count  ;

    // Flags for tracking the settings
    bool set_flux_flag = false ;
    bool set_e_range_flag = false ;
    bool set_e_center_flag = false ;
    bool set_n_count_flag = false ;
    bool divided_bin_flag = false ;

    std::vector<Zaki::Math::Quantity> t_bin_obs_set ;
    TH1F t_bin_hist ;
    size_t t_bin_chops = 0; 
};

//==============================================================
#endif /*DMSS_Bin_H*/
