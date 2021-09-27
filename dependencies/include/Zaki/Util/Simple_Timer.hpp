#ifndef Zaki_Util_Simple_Timer_H
#define Zaki_Util_Simple_Timer_H

#include <string>
#include <chrono>
// #include <algorithm>
#include <fstream>
#include <mutex>
#include <thread>

#include "Zaki/String/Directory.hpp"

//--------------------------------------------------------------
namespace Zaki::Util
{

//==============================================================
struct Operation
{
    std::string name ;
    double duration ;

    Operation(const std::string& in_name, double in_dur) 
        : name(in_name), duration(in_dur) {} 
};

//==============================================================
//                     Singleton design
class TimeManager
{

    //--------------------------------------------------------------
    private:

    std::string     m_sessionName   = "None";
    std::ofstream   m_OutputStream;
    int             m_ProfileCount = 0;
    std::mutex      m_lock;
    bool            m_activeSession = false;

    //--------------------------------------------------------------
    TimeManager() {}

    //--------------------------------------------------------------
    void IEndSession() ;

    //--------------------------------------------------------------
    void IBeginSession(const std::string&, const Zaki::String::Directory& = "Timer_Results.txt") ;

    //--------------------------------------------------------------
    public:

    //--------------------------------------------------------------
    ~TimeManager()
    {
        EndSession();
    }

    //--------------------------------------------------------------
    static void BeginSession(const std::string&, const Zaki::String::Directory& ="Timer_Results.txt") ;
    //--------------------------------------------------------------
    static void EndSession() ;

    //--------------------------------------------------------------
    void WriteProfile(const Operation& oper);

    //--------------------------------------------------------------
    void WriteHeader() ;

    //--------------------------------------------------------------
    void WriteFooter() ;

    //--------------------------------------------------------------
    static TimeManager& Get() ;
    //--------------------------------------------------------------
};

//==============================================================
class Timer
{
    //--------------------------------------------------------------
    private:

    Operation op ;
    std::chrono::time_point<std::chrono::high_resolution_clock> startTimePt ;
    bool stopped = false ;

    //--------------------------------------------------------------
    public:

    //--------------------------------------------------------------
    Timer(const std::string& in_scope) : op({in_scope, 0})
    {
        startTimePt = std::chrono::high_resolution_clock::now() ;
    }
    //--------------------------------------------------------------
    ~Timer()
    {
        if(!stopped) Stop() ; 
    }

    //--------------------------------------------------------------
    void Stop() ;
};

//==============================================================
//                  ostream << overloading
std::ostream& operator << (std::ostream &output, const Operation& o) ;

//--------------------------------------------------------------
} // End of namespace Zaki::Util
//==============================================================
//                    Profiling Macros
//==============================================================

#if defined(_MSC_VER)
#define Z_PRETTY_FUNCTION __FUNCSIG__
#else
#define Z_PRETTY_FUNCTION __PRETTY_FUNCTION__
#endif
//--------------------------------------------------------------
#define TIMING 1
#if TIMING
    #define Z_TIMER_SCOPE_HIDDEN(name, L) Zaki::Util::Timer __CONCAT(sim_timer,L)(name)
    #define Z_TIMER_SCOPE(name) Z_TIMER_SCOPE_HIDDEN(name, __COUNTER__)

/*
    // #define Z_TIMER_SCOPE(name) \
    // do { \
    //     Zaki::Util::Timer __CONCAT(sim_timer,__COUNTER__)(name); \
    //     return ; \
    // } while (0)
*/
    #define Z_TIMER() Z_TIMER_SCOPE(Z_PRETTY_FUNCTION)
#else
    #define Z_TIMER_SCOPE(name)
    #define Z_TIMER()
#endif

//==============================================================
#endif /*Zaki_Util_Simple_Timer_H*/
