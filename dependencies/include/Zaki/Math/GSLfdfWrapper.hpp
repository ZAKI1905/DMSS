#ifndef Zaki_Math_GSLfdfWrapper_H
#define Zaki_Math_GSLfdfWrapper_H

#include <gsl/gsl_math.h>

// Zaki::Util
#include "Zaki/Util/Logger.hpp"

//--------------------------------------------------------------
namespace Zaki::Math
{

  // double (* f) (double x, void * params);
  // double (* df) (double x, void * params);
  // void (* fdf) (double x, void * params, double * f, double * df);
  // void * params;

//==============================================================
template<typename FuncObj, typename MemFuncPtr, typename MemFdfPtr>
class GSLfdfWrapper : public gsl_function_fdf 
{
  public:
    /// Default Constructor
    GSLfdfWrapper()
    {
      f   = &GSLfdfWrapper::invoke_f ;
      df  = &GSLfdfWrapper::invoke_df ;
      fdf = &GSLfdfWrapper::invoke_fdf ;
      
      params   = this;
    }

    /// Constructor from Pointer-2-Obj and the member-function
    GSLfdfWrapper(FuncObj* objPtr, const MemFuncPtr& memFn, 
                  const MemFuncPtr& memDf, const MemFdfPtr& memFdf)
      : fObjPtr(objPtr), fMemFunc( memFn ), dfMemFunc( memDf ), 
        fdfMemFunc(memFdf)
    {
      f   = &GSLfdfWrapper::invoke_f ;
      df  = &GSLfdfWrapper::invoke_df ;
      fdf = &GSLfdfWrapper::invoke_fdf ;
      
      params   = this;
    }

    void SetMemberFunc(FuncObj* objPtr, const MemFuncPtr& memFn) 
    {
      fObjPtr   = objPtr ;
      fMemFunc  = memFn  ;
    }

  private:
    
    /// Pointer to the Object
    FuncObj* fObjPtr      = nullptr ;
    
    /// Function (f) method
    MemFuncPtr fMemFunc   = nullptr ;

    /// Derivative Function (df) method
    MemFuncPtr dfMemFunc  = nullptr ;

    /// Function & Derivative Function (fdf) method
    MemFdfPtr fdfMemFunc = nullptr ;

    double Eval_f(const double x) { return (fObjPtr->*fMemFunc)(x) ;}

    static double invoke_f(double x, void *params) 
    {
      return static_cast<GSLfdfWrapper*>(params)->Eval_f(x);
    }

    double Eval_df(const double x) { return (fObjPtr->*dfMemFunc)(x) ;}

    static double invoke_df(double x, void *params) 
    {
      return static_cast<GSLfdfWrapper*>(params)->Eval_df(x);
    }

    void Eval_fdf(const double x, double *f, double *df) 
    {
      (fObjPtr->*fdfMemFunc)(x, f, df) ;
    }

    static void invoke_fdf(double x, void *params, double *f, double *df) 
    {
      static_cast<GSLfdfWrapper*>(params)->Eval_fdf(x, f, df) ;
    }
};

//--------------------------------------------------------------
} // End of namespace Zaki::Math
//--------------------------------------------------------------
//==============================================================
#endif /*Zaki_Math_GSLfdfWrapper_H*/
