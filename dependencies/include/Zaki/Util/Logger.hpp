#ifndef Zaki_Util_Logger_H
#define Zaki_Util_Logger_H

#include <iostream>
#include <string>
#include <fstream>
#include <mutex>
#include <time.h>
#include <cstring>

#include "Zaki/String/String_Basic.hpp"
#include "Zaki/String/Directory.hpp"

//--------------------------------------------------------------
namespace Zaki::Util
{

// Forward declration
struct LogEntry ;

//==============================================================
enum class LogLevel
{
  Error = 0, Warning, Info, Verbose    
} ;

//==============================================================
// Singleton log manager class
class LogManager 
{

  public:

    // Destructor
    ~LogManager() ;

    // Copy constructor
    LogManager(const LogManager&) = delete ;

    static LogManager& Get();

    static void Emit(const LogEntry&);

    // Setters
    // static void SetSettings(const LogSettings&);
    static void SetLogLevels(const LogLevel&, const LogLevel&) ;
    static void SetLogLevels(const LogLevel&) ;
    static void SetLogFile(const Zaki::String::Directory&) ;
    static void SetBlackWhite(const bool&) ;

  private:
    // Constructor 1 with the default settings
    LogManager() { }

    void IEmit(const LogEntry&) ;
    void Print(const LogEntry&) ;
    void ShowReport() ;
    void AddFooter() ;
    void Export(const LogEntry&) ;
    std::string GetLogFileHeader() ;

    // Get current date/time, format is YYYY-MM-DD.HH:mm:ss
    std::string CurrentDateTime(const std::string&) const ;

    Zaki::String::Directory file_name = "";
    LogLevel file_thresh = LogLevel::Info ;
    LogLevel term_thresh = LogLevel::Info ;
    std::mutex mtx;
    
    bool export_flag = false ;
    bool first_time_writing = true ;

    // Black & white format
    bool blackWhite_format = false ;

    // LogSettings settings ; 
    int num_warnings = 0;
    int num_errors  = 0;
};

//==============================================================
struct LogFormat 
{
  LogLevel level ;
  std::string tag ;
  Zaki::String::Color tag_color, txt_color;

  // Default Constructor
  // LogFormat() ;

  // Constructor from LogLevel
  LogFormat(const LogLevel& in_lvl) 
    : level(in_lvl)
  {
    switch (in_lvl)
    {
    
    case LogLevel::Verbose:
        tag = "Note" ;
        tag_color = {Zaki::String::FGColor::Green} ;
        txt_color = {Zaki::String::FGColor::LGreen} ;
      break;
    
    case LogLevel::Info:
        tag = "Info" ;
        tag_color = {Zaki::String::FGColor::Cyan} ;
        txt_color = {Zaki::String::FGColor::LCyan} ;
      break;

    case LogLevel::Warning:
        tag = "Warning" ;
        tag_color = {Zaki::String::FGColor::Yellow} ;
        txt_color = {Zaki::String::FGColor::LYellow} ;
      break;

    case LogLevel::Error:
        tag = "Error" ;
        tag_color = {Zaki::String::FGColor::Red} ;
        txt_color = {Zaki::String::FGColor::LRed} ;
      break;

    default:
      break;
    }
  }

  static std::string ResetColor() 
  {
    return "\E[0m" ;
  }

};

//==============================================================
struct LogEntry
{
  friend class LogManager ;
  
  private:
    std::string Message ;
    int Line = 0 ;
    std::string FuncName ;
    std::string FileName ;
    // std::string time_stamp = "";
    LogFormat Format;


  public:
    // Constructor
    LogEntry(const std::string& in_mes, int in_line,
             const std::string& in_func, const std::string& in_file, const LogLevel& in_lev)
             : Message(in_mes), Line(in_line), FuncName(in_func),
               FileName(in_file), Format(in_lev)
    {
      LogManager::Emit(*this) ;
    }
} ;

//==============================================================

//--------------------------------------------------------------
} // End of namespace Zaki::Util
//==============================================================
//                     Logging
// #define __FILENAME__ (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)
#define __FILENAME__ strrchr("/" __FILE__, '/') + 1

#define Z_LOG_HIDDEN(MSG, LEVEL, L)  Zaki::Util::LogEntry __CONCAT(logEntry,L)(MSG,  __LINE__, __func__, __FILENAME__, LEVEL)
#define Z_LOG(MSG, LEVEL)  Z_LOG_HIDDEN(MSG, LEVEL, __COUNTER__)
/* #define Z_LOG(MSG, LEVEL) \
  // do { \
    // Zaki::Util::LogEntry __CONCAT(logEntry,__COUNTER__)(MSG,  __LINE__, __func__, __FILENAME__, LEVEL); \
    // return ; \
  // } while (0)
// #define Z_LOG(MSG, LEVEL)  Zaki::Util::LogEntry logEntry##__LINE__(MSG,  __LINE__, __func__, __FILE__, LEVEL)
 */

#define Z_LOG_NOTE(MSG)  Z_LOG(MSG, Zaki::Util::LogLevel::Verbose)
#define Z_LOG_INFO(MSG)  Z_LOG(MSG, Zaki::Util::LogLevel::Info)
#define Z_LOG_WARNING(MSG)  Z_LOG(MSG, Zaki::Util::LogLevel::Warning)
#define Z_LOG_ERROR(MSG)  Z_LOG(MSG, Zaki::Util::LogLevel::Error)

//==============================================================
#endif /*Zaki_Util_Logger_H*/
