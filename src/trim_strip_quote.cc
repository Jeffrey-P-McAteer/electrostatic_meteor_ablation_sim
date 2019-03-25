
#include <cstdio>
#include <string>

void trim_strip_quote(std::string &instr)
{
  
  size_t left,right;

  /* make sure it is long enough */
  if (instr.length() < 2) return;

  /* trim */
  left = instr.find_first_not_of(" \t\r\n\v\f");
  right = instr.find_last_not_of(" \t\r\n\v\f");
  instr = instr.substr(left,right-left+1);
  /* single quotes */
  if(instr.compare(0,1,"'") == 0 ) {
    if (instr.compare(instr.length()-1,1,"'") == 0) {
      instr.erase(0,1);
      instr.erase(instr.length()-1,1);
    } else {
      printf(";WARNING: error processing string, <%s>\n",instr.c_str());
    }
  }
  /* double quotes */
  else if (instr.compare(0,1,"\"") == 0 ) {
    if (instr.compare(instr.length()-1,1,"\"") == 0) {
      instr.erase(0,1);
      instr.erase(instr.length()-1,1);
    } else {
      printf(";WARNING: error processing string, <%s>\n",instr.c_str());
    }
  }
}
