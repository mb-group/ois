/* C3D: Computational protein Design by Duplication and Divergence (2021)
   Author: Duccio Malinverni, St.Jude Children's Research Hospital, Memphis, TN, USA 
   This file is covered by GPL-3.0 license (see LICENSE file in the root of this project.    
*/

#ifndef DEF_SELECTOR
#define DEF_SELECTOR

#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <algorithm>
#include "Utils.h"


void SEL_Select(std::string mutantsFile, float probThreshold, std::string outPrefix);
#endif
