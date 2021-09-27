#ifndef Zaki_Math_Func2D_H
#define Zaki_Math_Func2D_H

#include <memory>

//==============================================================
namespace Zaki::Math
{

//==============================================================
class GeneralObject
{
  public:
    GeneralObject() ;
    virtual GeneralObject* Clone() = 0 ;
    virtual ~GeneralObject() ;
};
//==============================================================
class Func2D
{
  protected:
    virtual Func2D* IClone() const = 0 ;
    // {
    //   return nullptr ;
    // };
    // virtual Func2D* IThis() = 0 ;


  public:
    Func2D() ;

    virtual ~Func2D() ;
    std::unique_ptr<Func2D> Clone() const 
    { return std::unique_ptr<Func2D>(IClone()); }

    // std::unique_ptr<Func2D> This() 
    // { return std::unique_ptr<Func2D>(IThis()); }

    virtual double Eval(const double x, const double y) = 0 ;
    std::string PtrStr() const ;
    // virtual void* GetObj() = 0 ;
    // virtual void UpdateObj(void*) = 0 ;

    // virtual GeneralObject* GetGenObj() = 0;

};
//==============================================================

//--------------------------------------------------------------
} // End of namespace Zaki::Math
//--------------------------------------------------------------

#endif /*Zaki_Math_Func2D_H*/