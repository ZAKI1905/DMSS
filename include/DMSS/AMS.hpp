#ifndef DMSS_AMS_H
#define DMSS_AMS_H

#include "DMSS/Satellite.hpp"

//==============================================================
class AMS : public Satellite
{
  public:
    
    /// Default constructor
    AMS() ;

    /// Copy constructor
    AMS(const AMS&) ;

    /// Destructor
    ~AMS() ;

    /// Acceptance (m^2*sr) of the detector as a function of energy (GeV)
    double Acceptance(double energy) const override ;

    /// Exposure time-fraction as a function of energy (GeV)
    double ExpTimeFrac(double energy) const override ;

  private:
    /// Overriden 'Satellite::IClone' method
    virtual AMS* IClone() const override ;
};



//==============================================================
#endif /*DMSS_AMS_H*/
