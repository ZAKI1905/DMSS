#ifndef Zaki_File_CSVIterator_H
#define Zaki_File_CSVIterator_H

#include "Zaki/File/CSVRow.hpp"

//--------------------------------------------------------------
namespace Zaki::File
{

//==============================================================
class CSVIterator
{

  //--------------------------------------------------------------
  public:

    // Constructors
    CSVIterator(std::istream& str, const char& in_delim=',') ;
    CSVIterator() ;

    void SetDelim(const char& in_delim);

    typedef std::input_iterator_tag     iterator_category;
    typedef CSVRow                      value_type;
    typedef std::size_t                 difference_type;
    typedef CSVRow*                     pointer;
    typedef CSVRow&                     reference;
    

    // Pre Increment
    CSVIterator& operator++() ;
    
    // Post increment
    CSVIterator operator++(int) ;  

    CSVRow const& operator*()   const ;     
    CSVRow const* operator->()  const ;   

    bool operator==(CSVIterator const& rhs) ;
    
    bool operator!=(CSVIterator const& rhs) ;

  //--------------------------------------------------------------

  private:
    std::istream*       m_str;
    CSVRow              m_row;

};

//--------------------------------------------------------------
} // End of namespace Zaki::File
//--------------------------------------------------------------

//==============================================================
#endif /*Zaki_File_CSVIterator_H*/
