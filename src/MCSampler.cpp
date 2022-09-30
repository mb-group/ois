/* OIS: Orthogonal Interacting Sequences (2022) 
   Author: Duccio Malinverni, St.Jude Children's Research Hospital, Memphis, TN, USA 
   This file is covered by GPL-3.0 license (see LICENSE file in the root of this project.    
*/

#include "MCSampler.h"

float MC_computeEnergy(float* __restrict__ P,uInt*  __restrict__ x,uInt  N,uInt q){
  float E=0;
  for(uInt i=0;i<N;i++){
    E+= P[h_at(i,x[i],q)];
    for(uInt j=i+1;j<N;j++)
      E+=P[J_at(i,j,x[i],x[j],N,q)];
  }
  return E;
}

void MC_sample(uInt* __restrict__ x, float* __restrict__ P,
	       uInt N,uInt q,uInt Nsweeps,uInt sweepsPerSample,
	       float* __restrict__ Es,uInt* __restrict__ data){
  
  // Setup MC variables
  uInt iFlip,qiNew;
  float deltaE;
  uInt nSaved = 0;
  
  // Perform the MC sampling
  float E=MC_computeEnergy(P,x,N,q);
  for(uInt sweep=0;sweep<Nsweeps;sweep++){
    for(uInt move=0;move<N;move++){
      iFlip = rand()%N;
      qiNew = rand()%q;

      // Compute the energy change
      deltaE = P[h_at(iFlip,qiNew,q)]-P[h_at(iFlip,x[iFlip],q)];
      for(uInt j=0;j<N;j++)
	if(iFlip!=j)
	  deltaE+=P[J_at(iFlip,j,qiNew,x[j],N,q)]-P[J_at(iFlip,j,x[iFlip],x[j],N,q)];
      
      // Accept with metropolis criterion
      if(deltaE<0 || std::exp(-deltaE)>(static_cast<float>(rand())/static_cast<float>(RAND_MAX))){
	x[iFlip] = qiNew;
	E += deltaE;
      }
    }

    // Record sample and energy
    if (sweep%sweepsPerSample == 0) {
      Es[nSaved] = E;
      for(uInt i=0;i<N;i++)
	data[nSaved*N+i]=x[i];
      nSaved++;
    }
  }
}
