#ifndef DMSS_HybPdf_H
#define DMSS_HybPdf_H

#include <gsl/gsl_integration.h>

// Root
#include <TF1.h>

#include "DMSS/Prog.hpp"

//==============================================================
class HybPdf : public Prog
{

  public:
    HybPdf();
    ~HybPdf() ;
    
    HybPdf(const HybPdf&) = delete ;
    HybPdf& operator=(const HybPdf&) = delete ;

    void SetPars(const std::vector<double>&) ;
    double IntegrandP(double);
    double IntegrandG(double);
    double NormalFactor(double x);
    double Integrate() ;
    double f(double, double) ;

  private:
    double rate, obs, obs_err;
    gsl_integration_workspace *GSL_workspace = nullptr ;
};

//==============================================================
#endif /*DMSS_HybPdf_H*/
