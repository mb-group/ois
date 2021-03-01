/* C3D: Computational protein Design by Duplication and Divergence (2021)
   Author: Duccio Malinverni, St.Jude Children's Research Hospital, Memphis, TN, USA 
   This file is covered by GPL-3.0 license (see LICENSE file in the root of this project.    
*/

#include "Utils.h"

std::vector<std::string> splitString(const std::string& s, const char& delimiter){
  std::string buff{""};
  std::vector<std::string> v;
  for(auto n:s)
    {
      if(n != delimiter) buff+=n; else
	if(n == delimiter && buff != "") { v.push_back(buff); buff = ""; }
    }
  if(buff != "") v.push_back(buff);
  
  return v;
}
