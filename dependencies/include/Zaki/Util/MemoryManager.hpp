#ifndef Zaki_Util_MemoryManager_H
#define Zaki_Util_MemoryManager_H

#include <iostream>
#include <string>
#include <fstream>
#include <mutex>

#include "Zaki/Util/Logger.hpp"
#include "Zaki/String/String_Basic.hpp"

//--------------------------------------------------------------
namespace Zaki::Util
{

// Forward declration
struct MemEntry ;

//==============================================================
// Singleton log manager class
class MemManager 
{

  public:

    // Destructor
    ~MemManager() ;

    // Copy constructor
    MemManager(const MemManager&) = delete ;

    static MemManager& Get();

    static void New(const MemEntry&);
    static void Delete(MemEntry&);


  private:

    // Constructor 1 with the default settings
    MemManager() { }

    void INew(const MemEntry&) ;
    void IDelete(MemEntry&) ;

    std::mutex mtx;
    static inline std::vector<MemEntry> heap_allocated = {} ; 
    long int total_mem_size = 0 ;

};

//==============================================================
struct MemEntry
{
  friend class MemManager ;

  public:
    enum class Mode
    {
      Delete=0, New
    };

  private:

    void Free()
    {
      Z_LOG_NOTE("Freeing heap allocated memory at: "+ Str() + ".");
      free(ptr);
    }

    void FailMessage()
    {
      Z_LOG_NOTE("Failed attempt to free a pointer at '" +
                  Str() + "' that wasn't allocated.");
    }

    std::string Str() const
    {
      std::stringstream ss;
      ss << ptr;
      return ss.str() ;
    }

    void* ptr ;
    size_t size = 0 ;
    // Mode mem_mode;

  public:
    // Constructor 1
    MemEntry(void* in_ptr, const Mode& in_mode)  
     : ptr(in_ptr)
    {
      if(in_mode == Mode::New)
        MemManager::New(*this) ;
      else
        MemManager::Delete(*this) ;
    }

    // Constructor 2
    MemEntry(void* in_ptr, const Mode& in_mode, const size_t& in_size)  
     : ptr(in_ptr), size(in_size)
    {
      if(in_mode == Mode::New)
        MemManager::New(*this) ;
      else
        MemManager::Delete(*this) ;
    }

    // MemEntry(const MemEntry&) = delete ;
    // MemEntry& operator=(const MemEntry&) = delete ;
    bool operator==(const MemEntry& other) const 
    {
      return ptr == other.ptr ;
    }


} ;

//==============================================================
} // End of namespace Zaki::Util
//==============================================================
//                     Memory Manager

#define Z_NEW_HIDDEN(PTR, SIZE, L)  Zaki::Util::MemEntry __CONCAT(memEntry,L)(PTR, Zaki::Util::MemEntry::Mode::New, SIZE)
#define Z_DELETE_HIDDEN(PTR, L)  Zaki::Util::MemEntry __CONCAT(memEntry,L)(PTR, Zaki::Util::MemEntry::Mode::Delete)

#define Z_NEW(PTR, SIZE)  Z_NEW_HIDDEN(PTR, SIZE, __COUNTER__)
#define Z_DELETE(PTR)  Z_DELETE_HIDDEN(PTR, __COUNTER__)

/* #define Z_NEW(PTR) \
//       do { \
//        Zaki::Util::MemEntry __CONCAT(memEntry,__COUNTER__)(PTR, Zaki::Util::MemEntry::Mode::New, SIZE); \
//        return ; \
//        } while (0)

// #define Z_DELETE(PTR) \
//       do { \
//        Zaki::Util::MemEntry __CONCAT(memEntry,__COUNTER__)(PTR, Zaki::Util::MemEntry::Mode::Delete); \
//        return ; \
//        } while (0)
 */
//==============================================================
#endif /*Zaki_Util_MemoryManager_H*/
