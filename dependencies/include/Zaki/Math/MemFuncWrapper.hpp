#ifndef MemFuncWrapper_H
#define MemFuncWrapper_H

#include "algorithm"

// Zaki::Math
#include "Zaki/Math/Func2D.hpp"

// Zaki::Util
#include "Zaki/Util/Logger.hpp"

// Zaki::Vector
#include "Zaki/Vector/Vector_Basic.hpp"

//--------------------------------------------------------------
namespace Zaki::Math
{

//==============================================================
template<typename FuncObj, typename MemFuncPtr >
class MemFuncWrapper : public Func2D
{

  public:
    MemFuncWrapper(const FuncObj& obj, const MemFuncPtr& memFn) 
      : fObj(obj), fMemFunc( memFn ) 
    {
      Z_LOG_NOTE("MemFuncWrapper constructor called from " + PtrStr() + ".") ;
    }
    
    MemFuncWrapper(const MemFuncWrapper &mfw2) 
    : fObj(mfw2.fObj), fMemFunc(mfw2.fMemFunc) 
    {
      Z_LOG_NOTE("MemFuncWrapper copy constructor: from " 
          + mfw2.PtrStr() + " --> " + PtrStr() + ".") ;

      // if constexpr(!std::is_copy_constructible_v<FuncObj>)
      //   Z_LOG_ERROR("Cannot copy 'MemFuncWrapper' because the wrapped object is not copy constructible.") ;

      cpy_const_called = true ;
    }

    ~MemFuncWrapper() 
    {
      Z_LOG_NOTE("MemFuncWrapper destructor called from " + PtrStr() + ".") ;
      // Z_LOG_NOTE(("heap_allocated: " +
      //             std::to_string(heap_allocated.size()) +
      //             ", heap_deleted: " + std::to_string(heap_deleted)
      //             + "."). ) ;
    }

    double Eval(const double x, const double y) override { return (fObj.*fMemFunc)( x , y) ;}

    // void* GetObj() override { return &fObj ;}
    FuncObj& GetObj() { return fObj ;}


    // GeneralObject* GetGenObj() override 
    // {
    //   return(dynamic_cast<GeneralObject*> (&fObj)) ;
    // }

    void UpdateObj(const FuncObj& otherfObj)
    {
      fObj = otherfObj;
    }

    // void UpdateObj(void* otherfObj) override
    // {
    //   fObj = *(FuncObj*)otherfObj;
    // }
    
    // Overloading CLass specific new operator
    static void* operator new(size_t sz)
    {
      void* m = malloc(sz);

      heap_allocated.push_back(m) ;
      return m;
    }
    // Overloading Class specific delete operator
    static void operator delete(void* m)
    {
      std::stringstream ss;
      ss << m;

      if (Zaki::Vector::Exists(m, heap_allocated)) 
      {
        Z_LOG_NOTE("Freeing heap allocated memory at: "+ ss.str() + ".");
        heap_allocated.erase(std::remove(heap_allocated.begin(), heap_allocated.end(), m)) ;
        free(m);
        // heap_deleted++ ;
      }
      else
      {
        Z_LOG_NOTE("Failed attempt to free a pointer at '" +
                    ss.str() + "' that wasn't allocated.");
      }
    }

  private:

    static inline std::vector<void*> heap_allocated = {} ; 
    // static inline unsigned int heap_deleted = 0; 

    virtual MemFuncWrapper* IClone() const override 
    {
      Z_LOG_NOTE("MemFuncWrapper::IClone() called from " + PtrStr() + ".") ;
      return (new MemFuncWrapper(*this)) ;
    }

    // virtual MemFuncWrapper* IThis() override
    // {
    //   return this ;
    // }

    bool cpy_const_called = false ;
    FuncObj fObj;
    MemFuncPtr fMemFunc ;

};

//--------------------------------------------------------------
} // End of namespace Zaki::Math
//--------------------------------------------------------------
//==============================================================
#endif /*MemFuncWrapper_H*/