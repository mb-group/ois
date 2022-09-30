/* OIS: Orthogonal Interacting Sequences (2022)  
   Author: Duccio Malinverni, St.Jude Children's Research Hospital, Memphis, TN, USA 
   This file is covered by GPL-3.0 license (see LICENSE file in the root of this project.    
*/

#ifndef DEF_EVALUATOR
#define DEF_EVALUATOR

#include "Utils.h"
#include "MCSampler.h"
# include <cmath>
#include <fstream>
#include <cstring>

void EVT_evaluate(Samples& samples, Parameters params, PottsModel model,float* &Es);
void EVT_computeDomainEnergy(float* __restrict__ P,uInt*  __restrict__ x,uInt  N,uInt q,
			      uInt nSplit,float &E_intra1,float &E_intra2, float &E_inter);
float EVT_computeDomainDeltaEnergy(float* __restrict__ P,uInt*  __restrict__ x,uInt  N,uInt q,uInt nSplit);
#endif
