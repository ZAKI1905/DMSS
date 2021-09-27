#ifndef Zaki_Math_GSLMultiFdfWrapper_H
#define Zaki_Math_GSLMultiFdfWrapper_H

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multiroots.h>

// Zaki::Util
#include "Zaki/Util/Logger.hpp"

//--------------------------------------------------------------
namespace Zaki::Math
{

//==============================================================
template< typename FuncObj, 
          typename MemFuncPtr, 
          typename MemDfPtr, 
          typename MemFdfPtr
        >
class GSLMultiFdfWrapper : public gsl_multiroot_function_fdf 
{
  public:
    /// Default Constructor
    GSLMultiFdfWrapper(const size_t& in_n)
    {
      f   = &GSLMultiFdfWrapper::invoke_f ;
      df  = &GSLMultiFdfWrapper::invoke_df ;
      fdf = &GSLMultiFdfWrapper::invoke_fdf ;

      n = in_n ;
      params   = this;

    }

    /// Constructor from Pointer-2-Obj and the member-function
    GSLMultiFdfWrapper(FuncObj* objPtr, const MemFuncPtr& memFn, 
                       const MemDfPtr& memDf, const MemFdfPtr& memFdf,
                       const size_t& in_n)
      : fObjPtr(objPtr), fMemFunc( memFn ), 
        dfMemFunc(memDf), fdfMemFunc(memFdf)
    {
      f   = &GSLMultiFdfWrapper::invoke_f ;
      df  = &GSLMultiFdfWrapper::invoke_df ;
      fdf = &GSLMultiFdfWrapper::invoke_fdf ;

      n = in_n ;
      params = this;
    }

    void SetMemberFunc(FuncObj* objPtr, const MemFuncPtr& memFn) 
    {
      fObjPtr   = objPtr ;
      fMemFunc  = memFn  ;
    }

  private:
    FuncObj* fObjPtr = nullptr ;

    /// Function (f) method
    MemFuncPtr fMemFunc   = nullptr ;

    /// Derivative Function (df) method
    MemDfPtr dfMemFunc  = nullptr ;

    /// Function & Derivative Function (fdf) method
    MemFdfPtr fdfMemFunc = nullptr ;

    int Eval_f(const gsl_vector* x, gsl_vector* f) 
    { 
      (fObjPtr->*fMemFunc)(x, f) ;

      return GSL_SUCCESS ;
    }

    static int invoke_f(const gsl_vector * x, void * p, gsl_vector * f) 
    {
      return static_cast<GSLMultiFdfWrapper*>(p)->Eval_f(x, f);
    }

    int Eval_df(const gsl_vector* x, gsl_matrix* J) 
    { 
      (fObjPtr->*dfMemFunc)(x, J) ;

      return GSL_SUCCESS ;
    }

    static int invoke_df(const gsl_vector * x, void * p, gsl_matrix * J) 
    {
      return static_cast<GSLMultiFdfWrapper*>(p)->Eval_df(x, J);
    }

    int Eval_fdf(const gsl_vector* x, gsl_vector*  f, gsl_matrix* J) 
    { 
      (fObjPtr->*fdfMemFunc)(x, f, J) ;

      return GSL_SUCCESS ;
    }

    static int invoke_fdf(const gsl_vector * x, void * p, gsl_vector* f, gsl_matrix * J) 
    {
      return static_cast<GSLMultiFdfWrapper*>(p)->Eval_fdf(x, f, J);
    }
};

//--------------------------------------------------------------
} // End of namespace Zaki::Math
//--------------------------------------------------------------
//==============================================================
#endif /*Zaki_Math_GSLMultiFdfWrapper_H*/
