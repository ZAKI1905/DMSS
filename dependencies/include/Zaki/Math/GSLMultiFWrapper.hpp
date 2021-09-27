#ifndef Zaki_Math_GSLMultiFWrapper_H
#define Zaki_Math_GSLMultiFWrapper_H

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>

// Zaki::Util
#include "Zaki/Util/Logger.hpp"

//--------------------------------------------------------------
namespace Zaki::Math
{

//==============================================================
template<typename FuncObj, typename MemFuncPtr >
class GSLMultiFWrapper : public gsl_multiroot_function 
{
  public:
    /// Default Constructor
    GSLMultiFWrapper(const size_t& in_n)
    {
      f = &GSLMultiFWrapper::invoke;
      n = in_n ;
      params   = this;
    }

    /// Constructor from Pointer-2-Obj and the member-function
    GSLMultiFWrapper(FuncObj* objPtr, const MemFuncPtr& memFn, const size_t& in_n)
      : fObjPtr(objPtr), fMemFunc( memFn )
    {
      f = &GSLMultiFWrapper::invoke;
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
    MemFuncPtr fMemFunc = nullptr ;

    int Eval(const gsl_vector* x, gsl_vector* f) 
    { 
      (fObjPtr->*fMemFunc)(x, f) ;

      return GSL_SUCCESS ;
    }

    static int invoke(const gsl_vector * x, void * p, gsl_vector * f) 
    {
      return static_cast<GSLMultiFWrapper*>(p)->Eval(x, f);
    }
};

//--------------------------------------------------------------
} // End of namespace Zaki::Math
//--------------------------------------------------------------
//==============================================================
#endif /*Zaki_Math_GSLMultiFWrapper_H*/
