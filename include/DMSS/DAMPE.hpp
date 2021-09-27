#ifndef DMSS_DAMPE_H
#define DMSS_DAMPE_H

#include "DMSS/Satellite.hpp"

//==============================================================
class DAMPE : public Satellite
{

  public:

    /// Default constructor
    DAMPE() ;

    /// Copy constructor
    DAMPE(const DAMPE&) ;

    /// Destructor
    ~DAMPE() ;

    /// Acceptance (m^2*sr) of the detector as a function of energy (GeV)
    double Acceptance(double energy) const override ;

  private:
    /// Overriden 'Satellite::IClone' method
    virtual DAMPE* IClone() const override ;
};



//==============================================================
#endif /*DMSS_DAMPE_H*/
