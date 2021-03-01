/* C3D: Computational protein Design by Duplication and Divergence (2021)
   Author: Duccio Malinverni, St.Jude Children's Research Hospital, Memphis, TN, USA 
   This file is covered by GPL-3.0 license (see LICENSE file in the root of this project.    
*/

#ifndef DEF_MCSAMPLER
#define DEF_MCSAMPLER

#include "Utils.h"
#include <cmath>

void MC_sample(uInt* __restrict__  x, float* __restrict__ P,
	       uInt N,uInt q,uInt Nsweeps,uInt sweepsPerSample,
	       float* __restrict__ Es,uInt* __restrict__ data);
float MC_computeEnergy(float* __restrict__ P,uInt* __restrict__ x,uInt N,uInt q);
#endif 
