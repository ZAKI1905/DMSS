#ifndef Zaki_Integrator_H
#define Zaki_Integrator_H

#include <gsl/gsl_integration.h>

#include "Zaki/Util/Logger.hpp"
#include "Zaki/Math/Math_Core.hpp"
#include "Zaki/Math/GSLFuncWrapper.hpp"

//--------------------------------------------------------------
namespace Zaki::Math
{

//==============================================================
class Integrator
{

  private:
    Quantity integral_result ;
    gsl_integration_workspace *GSL_workspace = nullptr ;
    size_t wrk_space_size ;
    double abs_err, rel_err ;

  public:

    Integrator(const size_t& in_wrk_space_size) 
    : wrk_space_size(in_wrk_space_size), 
      abs_err(1e-6), rel_err(1e-6)
    {
      GSL_workspace = gsl_integration_workspace_alloc(in_wrk_space_size) ;
    }

    Integrator(const size_t& in_wrk_space_size, 
    const double& in_eps) 
    : wrk_space_size(in_wrk_space_size), 
      abs_err(in_eps), rel_err(in_eps)
    {
      GSL_workspace = gsl_integration_workspace_alloc(in_wrk_space_size) ;
    }

    Integrator(const size_t& in_wrk_space_size, 
    const double& in_epsabs,
    const double& in_epsrel ) 
    : wrk_space_size(in_wrk_space_size), 
      abs_err(in_epsabs), rel_err(in_epsrel)
    {
      GSL_workspace = gsl_integration_workspace_alloc(in_wrk_space_size) ;
    }

    virtual ~Integrator() 
    {
      gsl_integration_workspace_free(GSL_workspace);
    }
    
    Integrator(const Integrator&) = delete ;
    Integrator& operator=(const Integrator&) = delete ;

    virtual double Integrand(double) const = 0 ;
    virtual Integrator* GetIntegratorPtr() = 0 ;

    Quantity GetIntegrationResult() const
    {
      return integral_result ;
    }

    Quantity Integrate(double min, double max) 
    {
      GSLFuncWrapper<Integrator, double (Integrator::*)(double) const> 
          tmp_f(GetIntegratorPtr(), &Integrator::Integrand);     

      gsl_function tmp_f_GSL = static_cast<gsl_function> (tmp_f) ; 

      gsl_integration_qag(  &tmp_f_GSL, min, max, abs_err, rel_err, 
                            wrk_space_size, 1,
                            GSL_workspace, &integral_result.val,
                            &integral_result.err);

      return integral_result ;
    }
};

//--------------------------------------------------------------
} // End of namespace Zaki::Math
//--------------------------------------------------------------
//==============================================================
#endif /*Zaki_Integrator_H*/
