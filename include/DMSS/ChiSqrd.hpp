#ifndef DMSS_ChiSqrd_H
#define DMSS_ChiSqrd_H

#include <Minuit2/FCNBase.h>

#include <Zaki/File/CSVIterator.hpp>

#include "DMSS/Bin.hpp"

//==============================================================
namespace ROOT {

namespace Minuit2 {


class ChiSqrd : public FCNBase, public Prog
{

  public:
    
    /// Constructor
    ChiSqrd() : Prog("ChiSqrd", true) {}
    
    /// Copy constructor
    ChiSqrd(const ChiSqrd& other)
    : Prog("ChiSqrd", true),
    theErrorDef(other.theErrorDef),
    obs_set(other.obs_set)
    {
        for (auto b: sig_shape_hist) {
            b = *(static_cast<TH1F*>(b.Clone())) ;
        }
    }
    
    virtual double operator()(const std::vector<double>& par) const override ;
    
    ~ChiSqrd() {}

    double Up() const override {return theErrorDef;}

    void SetErrorDef(double def) override
    { 
      theErrorDef = def; 
      char tmp[100] ;
      sprintf(tmp, "%.3f", def) ;
      Z_LOG_INFO("The error definition is set to: " + std::string(tmp) + ".") ;
    }

    void AddObsCounts(const Bin&) ;
    void AddSigShape(const TH1F&) ;

  private:
    double theErrorDef = 1 ;
    std::vector<std::vector<Zaki::Math::Quantity> > obs_set ;
    std::vector<TH1F> sig_shape_hist ;
};

  }  // namespace Minuit2

}  // namespace ROOT
//==============================================================
#endif /*DMSS_ChiSqrd_H*/
