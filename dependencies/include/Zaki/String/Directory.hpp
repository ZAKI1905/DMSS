#ifndef Zaki_String_Directory_H
#define Zaki_String_Directory_H

#include "string"
//--------------------------------------------------------------
namespace Zaki::String
{

//==============================================================

class Directory
{

  public:
    Directory(const std::string& in_path) 
      : full_path(in_path)
    {
      // if(*(full_path.end()-1) == '/')
      // {
      //   // strip the extra backslash from the end
      //   full_path = full_path.substr(0, full_path.find_last_of("/")) ;
      // }

      // if(*full_path.begin() == '/')
      // {
      //   // strip the extra backslash from the beginning
      //   full_path = full_path.substr(1) ;
      // }
    };

    // Overloading for const char* input
    Directory(const char* in_path) 
      : full_path(in_path)
    {
      // if(*(full_path.end()-1) == '/')
      //   // strip the extra backslash from the end
      //   full_path = full_path.substr(0, full_path.find_last_of("/")) ;

      // if(*full_path.begin() == '/')
      // {
      //   // strip the extra backslash from the beginning:
      //   full_path = full_path.substr(1) ;
      // }
    };

    // Setters
    void SetPath(const std::string&) ;

    // Getters
    static Directory ThisFileDir(const char*) ;
    static Directory ThisFile(const char*) ;
    Directory NoExt() const ;
    Directory ParentDir() const ;
    std::string GetPath() const ;
    std::string Str() const ;
    void Print() const ;

    Directory operator+(const std::string&) const;
    Directory operator+(const char*) const;
    Directory operator+(const Directory&) const;
    // Overloading postfix -- operator
    Directory operator--(int) ;

  private:
    std::string full_path = "";

};
//------------------------------------------------------------
//==============================================================
//                  ostream << overloading
std::ostream& operator << (std::ostream &output, const Zaki::String::Directory&) ;

//==============================================================
} // End of namespace Zaki::String
//--------------------------------------------------------------

//==============================================================
#endif /*Zaki_String_Directory_H*/