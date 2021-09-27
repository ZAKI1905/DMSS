#ifndef Zaki_String_TextBox_H
#define Zaki_String_TextBox_H

#include "sstream"
#include "vector"

#include "Zaki/String/String_Basic.hpp"

//--------------------------------------------------------------
namespace Zaki::String
{

//==============================================================
class TextBox
  {
  
  public:
    // Constructor 1
    TextBox();
    // Constructor 2
    TextBox(const std::vector<std::string>&);

    enum TextAlignment
    {
      left = 0, center = 1, right = 2
    };

    void SetText(const std::vector<std::string>&) ;
    void AddText(const std::string&) ;
    void AddText(const std::vector<std::string>&) ;

    void SetAlignment(const TextAlignment&) ;
    void SetFrameColor(const Color&) ;
    void SetTextColor(const Color&) ;
    void SetPadColor(const Color&) ;
    void SetPadLeft(const size_t) ;
    void SetPadRight(const size_t) ;
    void SetPadding(const size_t, const int = -1) ;
    
    // Clears the screen first & moves the cursor to (1,1)
    void EnableClearScreen() ;

    std::string Str() ;
    void Print() ;
    void MakeFrame() ;
    void MakeStrFrame() ;

  private:
    const char* top_left_corn = "\u2554";
    const char* top_right_corn = "\u2557";
    const char* bottom_left_corn = "\u255A";
    const char* bottom_right_corn = "\u255D";
    const char* horizontal_line = "\u2550";
    const char* vertical_line = "\u2551";

    Color reset_color; // "\E[0m"

    Color frame_color;
    Color text_color;
    Color pad_color;

    // Clears the screen first & moves the cursor to (1,1)
    bool clear_screen = false; 

    size_t height = 0, width = 0, l_pad = 1, r_pad = 1 ;

    TextAlignment text_align = left ;

    std::stringstream top_frame, mid_section, bottom_frame;
    std::vector<std::string> text ;

    // std::string Multiply(const char*, size_t) ;
    size_t CorrectEncodedChars(const std::string&) ;

  };
  //------------------------------------------------------------

//--------------------------------------------------------------
} // End of namespace Zaki
//--------------------------------------------------------------

//==============================================================
#endif /*Zaki_String_TextBox_H*/