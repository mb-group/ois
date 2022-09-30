/* OIS: Orthogonal Interacting Sequences (2022)  
   Author: Duccio Malinverni, St.Jude Children's Research Hospital, Memphis, TN, USA 
   This file is covered by GPL-3.0 license (see LICENSE file in the root of this project.    
*/

#ifndef DEF_ORTHOGENERATOR
#define DEF_ORTHOGENERATOR

#include "Utils.h"
#include "Evaluator.h"
#include <algorithm>
#include <fstream>
#include <cmath>
#include <cstring>
#include <iomanip> 

void ORT_generate(Samples native, Parameters params, PottsModel model);
#endif
