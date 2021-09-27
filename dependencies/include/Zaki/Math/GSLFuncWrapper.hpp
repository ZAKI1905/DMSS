#ifndef Zaki_Math_GSLFuncWrapper_H
#define Zaki_Math_GSLFuncWrapper_H

#include <gsl/gsl_math.h>

// Zaki::Util
#include "Zaki/Util/Logger.hpp"

//--------------------------------------------------------------
namespace Zaki::Math
{

//==============================================================
template<typename FuncObj, typename MemFuncPtr >
class GSLFuncWrapper : public gsl_function 
{
  public:
    /// Default Constructor
    GSLFuncWrapper()
    {
      function = &GSLFuncWrapper::invoke;
      params=this;
    }

    /// Constructor from Pointer-2-Obj and the member-function
    GSLFuncWrapper(FuncObj* objPtr, const MemFuncPtr& memFn)
      : fObjPtr(objPtr), fMemFunc( memFn )
    {
      function = &GSLFuncWrapper::invoke;
      params=this;
    }

    void SetMemberFunc(FuncObj* objPtr, const MemFuncPtr& memFn) 
    {
      fObjPtr   = objPtr ;
      fMemFunc  = memFn  ;
    }

  private:
    FuncObj* fObjPtr = nullptr ;
    MemFuncPtr fMemFunc = nullptr ;
    double Eval(const double x) { return (fObjPtr->*fMemFunc)(x) ;}

    static double invoke(double x, void *params) {
      return static_cast<GSLFuncWrapper*>(params)->Eval(x);
  }
};

//--------------------------------------------------------------
} // End of namespace Zaki::Math
//--------------------------------------------------------------
//==============================================================
#endif /*Zaki_Math_GSLFuncWrapper_H*/
