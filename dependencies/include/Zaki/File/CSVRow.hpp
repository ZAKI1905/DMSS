#ifndef Zaki_File_CSVRow_H
#define Zaki_File_CSVRow_H

#include <iterator>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

//--------------------------------------------------------------
namespace Zaki::File
{

//==============================================================
class CSVRow
{

  //--------------------------------------------------------------
  public:

    CSVRow() ;

    std::string const& operator[](std::size_t index) const ;

    std::size_t size() const ;
      
    void readNextRow(std::istream& str) ;
    void SetDelim(const char& in_delim) ;

  //--------------------------------------------------------------
  private:
    std::vector<std::string>    m_data;
    char delim = ',' ;
};

//--------------------------------------------------------------
} // End of namespace Zaki::File
//--------------------------------------------------------------

std::istream& operator>>(std::istream& str, Zaki::File::CSVRow& data) ;
//==============================================================
#endif /*Zaki_File_CSVRow_H*/
