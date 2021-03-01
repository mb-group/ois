/* C3D: Computational protein Design by Duplication and Divergence (2021)
   Author: Duccio Malinverni, St.Jude Children's Research Hospital, Memphis, TN, USA 
   This file is covered by GPL-3.0 license (see LICENSE file in the root of this project.    
*/

#ifndef DEF_TRAINING
#define DEF_TRAINING

#include "MCSampler.h"
#include <mpi.h>

void TRN_train(float* P,Samples samples,Parameters params);
void TRN_computeFrequencies(float* F, uInt* data, uInt Nsamples, uInt N, uInt q);

#endif
