#ifndef Zaki_Util_Profile_Timer_H
#define Zaki_Util_Profile_Timer_H

#include <string>
#include <chrono>
#include <vector>
#include <fstream>
#include <mutex>
#include <thread>

#include "Zaki/String/Directory.hpp"

//--------------------------------------------------------------
namespace Zaki::Util
{
//==============================================================
struct Profile
{
    std::string name ;
    double duration ;
    std::vector<int> idx ;

    // Constructor 1
    Profile(const std::string& in_name, double in_dur,
            const std::vector<int> in_idx) 
        : name(in_name), duration(in_dur), idx(in_idx) {} 
};

//==============================================================
//                     Singleton design
class TimeProfileManager
{

    //--------------------------------------------------------------
    private:

    std::string     m_SessionName   = "None";
    std::ofstream   m_OutputStream;
    int             m_ProfileCount = 0;
    std::mutex      m_Lock;
    bool            m_ActiveSession = false;
    bool            m_ActiveSpecialSession = false ;
    std::vector<double> special_values ;

    //--------------------------------------------------------------
    TimeProfileManager() {}
    //--------------------------------------------------------------
    void IEndSession() ;
    //--------------------------------------------------------------
    void IBeginSession(const std::string&, const Zaki::String::Directory& = "Time_Profile_Results.dat") ;
    //--------------------------------------------------------------
    void IEndSpecialSession() ;
    //--------------------------------------------------------------
    void IBeginSpecialSession(const std::string&, const Zaki::String::Directory& = "Time_Profile_Special_Results.dat") ;
    //--------------------------------------------------------------
    std::vector<double> IGetSpecialValues() ;
    //--------------------------------------------------------------
    
    public:

    //--------------------------------------------------------------
    ~TimeProfileManager()
    {
        EndSession();
    }

    //--------------------------------------------------------------
    static void BeginSession(const std::string&, const Zaki::String::Directory& = "Time_Profile_Results.dat") ;
    //--------------------------------------------------------------
    static void EndSession() ;
    //--------------------------------------------------------------
    static void BeginSpecialSession(const std::string&, const Zaki::String::Directory& ="Time_Profile_Special_Results.dat") ;
    //--------------------------------------------------------------
    static void EndSpecialSession() ;
    //--------------------------------------------------------------
    static std::vector<double> GetSpecialValues() ;
    //--------------------------------------------------------------
    void WriteProfile(const Profile& oper);
    //--------------------------------------------------------------
    void WriteHeader() ;
    //--------------------------------------------------------------
    void WriteFooter() ;
    //--------------------------------------------------------------
    static TimeProfileManager& Get() ;
    //--------------------------------------------------------------
};

//==============================================================
class ProfileTimer
{
    //--------------------------------------------------------------
    private:

    Profile prof ;
    std::chrono::time_point<std::chrono::high_resolution_clock> startTimePt ;
    bool stopped = false ;

    //--------------------------------------------------------------
    public:

    //--------------------------------------------------------------
    ProfileTimer(const std::string& in_scope, const std::vector<int>& in_idx) 
        : prof(in_scope, 0, in_idx)
    {
        startTimePt = std::chrono::high_resolution_clock::now() ;
    }
    
    //--------------------------------------------------------------
    ProfileTimer(const std::string& in_scope) 
        : prof(in_scope, 0, {})
    {
        startTimePt = std::chrono::high_resolution_clock::now() ;
    }
    
    //--------------------------------------------------------------
    ~ProfileTimer() ;

    //--------------------------------------------------------------
    void Stop() ;
};

//==============================================================
//                  ostream << overloading
std::ostream& operator << (std::ostream &output, Profile p) ;

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
#define PROF_TIMING 1
#if PROF_TIMING
    #define Z_TIME_PROFILE_HIDDEN(name, set, L) Zaki::Util::ProfileTimer __CONCAT(prof_timer,L)(name, set)
    #define Z_TIME_PROFILE_SIMPLE_HIDDEN(name, L) Zaki::Util::ProfileTimer __CONCAT(prof_timer,L)(name)
    #define Z_TIME_PROFILE(name, set) Z_TIME_PROFILE_HIDDEN(name, set, __COUNTER__)
    #define Z_TIME_PROFILE_SIMPLE(name) Z_TIME_PROFILE_SIMPLE_HIDDEN(name, __COUNTER__) 
/*
    // #define Z_TIME_PROFILE(name, set) \
    // do { \
    //     Zaki::Util::ProfileTimer __CONCAT(prof_timer, __COUNTER__)(name, set); \
    //     return ; \
    // } while (0)
    
    // #define Z_TIME_PROFILE_SIMPLE(name) \
    // do { \
    //     Zaki::Util::ProfileTimer Zaki::Util::ProfileTimer __CONCAT(prof_timer,__COUNTER__)(name); \
    //     return ; \
    // } while (0)
*/
#else
    #define Z_TIME_PROFILE(name, set)
    #define Z_TIME_PROFILE_SIMPLE(name)
#endif

//==============================================================
#endif /*Zaki_Util_Profile_Timer_H*/
