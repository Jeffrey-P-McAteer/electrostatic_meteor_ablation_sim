// An enhanced map class, mapping string to string. New features include easy
// assignment and conversion to and from string values.

// Written by Yann Tambouret starting in 2010


#ifndef DICTIONARY_H
#define DICTIONARY_H

#include <string>
#include <map>
#include <sstream>

//using namespace std;

class dictionary : public std::map< std::string, std::string>
{
  
 protected:
  
 public:
  
  //  Constructors:
  //  -- inherit from map class

  // Copy Constructor
  //  -- inherit from map class

  //  Assignment
  const dictionary& assign(const char* key,const char* val){
    std::ostringstream tmpstring;
    tmpstring << val;
    this->insert(std::pair<std::string,std::string>(key,tmpstring.str()));
    return *this;
  }

  const dictionary& assign(const char* key,int val){
    std::ostringstream tmpstring;
    tmpstring << val;
    this->insert(std::pair<std::string,std::string>(key,tmpstring.str()));
    return *this;
  }

  const dictionary& assign(const char* key,float val){
    std::ostringstream tmpstring;
    tmpstring << val;
    this->insert(std::pair<std::string,std::string>(key,tmpstring.str()));
    return *this;
  }

  const dictionary& assign(const char* key,double val){
    std::ostringstream tmpstring;
    tmpstring << val;
    this->insert(std::pair<std::string,std::string>(key,tmpstring.str()));
    return *this;
  }
  //  Destructor
  //  -- inherit from map class
  
  //  Access and setting -- Speed is everything

  inline int toint(const char* key)  {
    return atoi(this->find(key)->second.c_str());
  }

  inline float toflt(const char* key)  {
    return atof(this->find(key)->second.c_str());
  }

  inline double todbl(const char* key)  {
    std::istringstream iss(this->find(key)->second);
    double tmp=0;
    if (!(iss >> std::dec >> tmp).fail())
      return tmp;
  }

};

#endif // ArrayNd_ranged_H 
