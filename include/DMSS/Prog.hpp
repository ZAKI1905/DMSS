#ifndef DMSS_Prog_H
#define DMSS_Prog_H

#include "DMSS/DMSSConfig.h"

#include <Zaki/Util/Logger.hpp>
#include <Zaki/Util/Instrumentor.hpp>
#include <Zaki/Util/Simple_Timer.hpp>
#include <Zaki/String/Directory.hpp>

#include "stdio.h"
#include <atomic>

#define DMSS_PROG_DEBUG_MODE 1

//==============================================================
class Prog
{
  //--------------------------------------------------------------
  protected:

    /// Prints the banned of the program
    void ShowBanner();

    /// The number of derived object instances
    static inline std::atomic<size_t> counter = 0;

    /// The work directory
    Zaki::String::Directory wrk_dir = "" ;

    /// The label of the class
    std::string name        = "" ;

    // Flags

    /// Returns if the name is set
    bool set_name_flag = false    ;

    /// Returns if the work directory is set
    bool set_wrk_dir_flag = false ;
    
  //--------------------------------------------------------------
  private:

#if DMSS_PROG_DEBUG_MODE
    /// Track object creation and destruction
    bool track_objs_flag = false ;
#endif
  //--------------------------------------------------------------
  public:

    /// Default Constructor for faster object creation
    ///  -> No name assignment or object tracking
    Prog() ;
    
    /// Constructor via a name and no tracking options
    /// -> Name is set but no object tracking
    Prog(const std::string&) ;

    /// Constructor via a name and tracking options
    /// This is slower and shouldn't be used for objects 
    /// that are created thousands of time, i.e. particles, etc.
    /// -> Name is set and object tracking is on
    Prog(const std::string&, const bool&) ;

    /// Destructor
    virtual ~Prog() ;

    /// Overloading CLass specific new operator
    static void* operator new(size_t sz) ;

    /// Overloading CLass specific delete operator
    static void operator delete(void* m);

    // Setters

    /// Sets the work directory
    virtual void SetWrkDir(const Zaki::String::Directory&) ;

    /// Sets the name of the class
    virtual void SetName(const std::string&) ;

    // void SetObjTracking(const bool) ;

    // Getters

    /// Returns the work directory
    virtual Zaki::String::Directory GetWrkDir() const ;

    /// Returns the name of the class
    virtual std::string GetName() const ;

    /// Print method for the derived classes
    virtual void Print() const ;

    /// String form of the pointer to the object
    std::string PtrStr() const ;

};


//==============================================================
#endif /*DMSS_Prog_H*/
