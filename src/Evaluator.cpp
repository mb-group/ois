#include "Evaluator.h"

void EVT_evaluate(Samples& samples, Parameters params,PottsModel model,float* &Es){
  if (params.mutantsListFile==""){
    // Evaluate energies on all provided samples   
    Es = new float[4*samples.Nsamples];
    float E_intra1,E_intra2,E_inter;
    for (uInt i=0;i<samples.Nsamples;i++){
      EVT_computeDomainEnergy(model.P,&(samples.data[samples.N*i]),samples.N,samples.q,
			      params.nDomainSplit,E_intra1,E_intra2,E_inter);
      Es[4*i]=E_intra1+E_intra2+E_inter;
      Es[4*i+1]=E_intra1;
      Es[4*i+2]=E_intra2;
      Es[4*i+3]=E_inter;
    }
  }
  else if(params.mutantsListFile!=""){
    // Compute all mutants provided by a mutation list
    std::vector<float> tmpEs;
    std::ifstream in(params.mutantsListFile.c_str());
    std::string tmp;
    std::vector<std::string> splitLine,splitMut;
    uInt* mutant=new uInt[model.N];
    float E_intra1,E_intra2,E_inter;
    
    while(!in.eof()){
      getline(in,tmp);
      if(tmp!=""){
	splitLine=splitString(tmp,' ');
	std::memcpy(mutant,&samples.data[0],model.N*sizeof(uInt));
	for(unsigned int i=0;i<splitLine.size();i++){
	  splitMut=splitString(splitLine[i],'_');
	  if(mutant[atoi(splitMut[1].c_str())]!=atoi(splitMut[0].c_str())){
	    std::cout<<"Inconsistent mutation, native symbol different."<<std::endl;
	    std::cout<<splitLine[i]<<std::endl;
	    exit(1);
	  }
	  mutant[atoi(splitMut[1].c_str())]=atoi(splitMut[2].c_str());
	}
	EVT_computeDomainEnergy(model.P,mutant,model.N,model.q,
				params.nDomainSplit,E_intra1,E_intra2,E_inter);
	tmpEs.push_back(E_intra1+E_intra2+E_inter);
	tmpEs.push_back(E_intra1);
	tmpEs.push_back(E_intra2);
	tmpEs.push_back(E_inter);
      }
    }
    Es= new float[tmpEs.size()];
    std::memcpy(Es,&tmpEs[0],tmpEs.size()*sizeof(float));
    samples.Nsamples=tmpEs.size()/4;
  }
}


void EVT_computeDomainEnergy(float* __restrict__ P,uInt*  __restrict__ x,uInt  N,uInt q,
			      uInt nSplit,float &E_intra1,float &E_intra2, float &E_inter){
  E_intra1=0;
  E_intra2=0;
  E_inter=0;
  
  for(uInt i=0;i<=nSplit;i++){
    E_intra1+= P[h_at(i,x[i],q)];
    for(uInt j=i+1;j<=nSplit;j++)
      E_intra1+=P[J_at(i,j,x[i],x[j],N,q)];
    
    for(uInt j=nSplit+1;j<N;j++)
      E_inter+=P[J_at(i,j,x[i],x[j],N,q)];
  }
  
  for(uInt i=nSplit+1;i<N;i++){
    E_intra2+= P[h_at(i,x[i],q)];
    for(uInt j=i+1;j<N;j++)
      E_intra2+=P[J_at(i,j,x[i],x[j],N,q)];
  }
}

float EVT_computeDomainDeltaEnergy(float* __restrict__ P,uInt*  __restrict__ x,uInt  N,uInt q){
  return 0;
}
