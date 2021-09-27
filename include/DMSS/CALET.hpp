#ifndef DMSS_CALET_H
#define DMSS_CALET_H

#include "DMSS/Satellite.hpp"

//==============================================================
class CALET : public Satellite
{

  public:

    /// Default constructor
    CALET() ;

    /// Copy constructor
    CALET(const CALET&) ;

    /// Destructor
    ~CALET() ;

    /// Acceptance (m^2*sr) of the detector as a function of energy (GeV)
    double Acceptance(double energy) const override ;

    /// Exposure time-fraction as a function of energy (GeV)
    double ExpTimeFrac(double energy) const override ;

  private:
    /// Overriden 'Satellite::IClone' method
    virtual CALET* IClone() const override ;
};



//==============================================================
#endif /*DMSS_CALET_H*/
