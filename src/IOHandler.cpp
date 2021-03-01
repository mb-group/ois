/* C3D: Computational protein Design by Duplication and Divergence (2021)
   Author: Duccio Malinverni, St.Jude Children's Research Hospital, Memphis, TN, USA 
   This file is covered by GPL-3.0 license (see LICENSE file in the root of this project.    
*/

#include "IOHandler.h"
using namespace std;

IOHandler::IOHandler()
{}

IOHandler::~IOHandler()
{}
void IOHandler::setOutPrefix(std::string out){
  outPrefix = out;
}

Parameters IOHandler::loadArgsTrain(int argc, char** argv){
  //If help option or nothing passed as arguments, display help and exit
  if(argc==2 || !string(argv[2]).compare("-h") || !string(argv[2]).compare("--help")){
    cout<<"Usage: c3d train -f samplesFile [options]"<<endl;
    cout<<"          -f       : Input file for learning (space delimited raw format)"<<endl;
    cout<<"          -N       : Total number of sweeps to perform [default 10000]"<<endl;
    cout<<"          -n       : Number of sweeps between recording two samples [default 10]"<<endl;
    cout<<"          -p       : Potts model starting parameters file in prm format"<<endl;
    cout<<"          -o       : Output prefix for saving files [default \"output\"]"<<endl;
    cout<<"          -b       : Number of iterations of the BLM optimization [default 1000]"<<endl;
    cout<<"          -l       : Regularization parameter lambda [default 0.01]"<<endl;
    cout<<"          -r       : Learning rate eta [default 0.01]"<<endl;
    cout<<"          --seed   : Random number generator seed [default 0 == time]."<<endl;
    exit(0);
  }

  Parameters params=parseInputArgs(argc,argv);
  return params;
}

Parameters IOHandler::loadArgsOrtho(int argc, char** argv){
  //If help option or nothing passed as arguments, display help and exit
  if(argc==2 || !string(argv[2]).compare("-h") || !string(argv[2]).compare("--help")){
    cout<<"Usage: c3d generate -f nativeFile -p prmFile [options]"<<endl;
    cout<<"          -f       : Native sample file (space delimited one-line sample file)"<<endl;
    cout<<"          -p       : Potts model parameters file in prm format"<<endl;
    cout<<"          -M       : Number of mutants to compute. [default 100]"<<endl;
    cout<<"          -m       : Number of point mutations per mutant [default 0]"<<endl;
    cout<<"          -m2      : Number of point mutations per mutant on the second domain [default 0]"<<endl;
    cout<<"          -T       : Virtual scaling temperature in the MCMC sampling [default 1]."<<endl;
    cout<<"          -o       : Output prefix for saving files [default \"output\"]"<<endl;
    cout<<"          -pi      : List of positions (0-based indexes) onto which to restrict the mutations [default None]"<<endl;
    cout<<"          -ns      : Last index (inclusive, 0-based) of first domain, \n                  "
                                "used for domain split E computation.  [default None]"<<endl;
    cout<<"          -n       : Number of sweeps between recording two mutants [default 10]"<<endl;
    cout<<"          --seed   : Random number generator seed [default 0 == time]."<<endl;
    exit(0);
  }

  Parameters params=parseInputArgs(argc,argv);
  return params;
}

Parameters IOHandler::loadArgsSelect(int argc, char** argv){
  //If help option or nothing passed as arguments, display help and exit
  if(argc==2 || !string(argv[2]).compare("-h") || !string(argv[2]).compare("--help")){
    cout<<"Usage: c3d select -f mutantsFile -t probThreshold [options]]"<<endl;
    cout<<"          -f       : Mutants file, as output by c3d generate (comprising mutations and energies)."<<endl;
    cout<<"          -t       : Probability threshold (in [0,1]) to select orthogonal mutants."<<endl;
    cout<<"          -o       : Output prefix for saving files [default \"output\"]"<<endl;
    exit(0);
  }

  Parameters params=parseInputArgs(argc,argv);
  return params;
}


Parameters IOHandler::parseInputArgs(int argc, char** argv){
  //Define default parameters
  Parameters params;
  params.randomSeed=0;
  params.Nsweeps = 10000;
  params.sweepsPerSample=10;
  params.outPrefix = "output";  
  params.prmInputFile=""; 

  params.samplesInputFile="";
  params.blmIterations=1000;
  params.regularizationLambda = 0.01;
  params.learningRate = 0.01;

  params.nPointMutations=0;
  params.nPointMutationsDom2=0;
  params.nMutants=100;
  params.T=1.;
  params.probThreshold=-1;
    
  params.mutantsListFile="";
  params.nDomainSplit=0;

  params.mutatePositionsFile="";
  // Parse command line arguments
  for (uInt i=0;i<argc;i++){
    if(!string(argv[i]).compare("-N")){
      params.Nsweeps = atoi(argv[i+1]);
      i++;
    }
    if(!string(argv[i]).compare("-n")){
      params.sweepsPerSample = atoi(argv[i+1]);
      i++;
    }
    if(!string(argv[i]).compare("-f")){
      params.samplesInputFile = string(argv[i+1]);
      i++;
    }
    if(!string(argv[i]).compare("-p")){
      params.prmInputFile = string(argv[i+1]);
      i++;
    }
    if(!string(argv[i]).compare("-o")){
      params.outPrefix = string(argv[i+1]);
      i++;
    }
    if(!string(argv[i]).compare("-b")){
      params.blmIterations = atoi(argv[i+1]);
      i++;
    }
    if(!string(argv[i]).compare("-l")){
      params.regularizationLambda = atof(argv[i+1]);
      i++;
    }
    if(!string(argv[i]).compare("-r")){
      params.learningRate = atof(argv[i+1]);
      i++;
    }
    if(!string(argv[i]).compare("--seed")){
      params.randomSeed = atoi(argv[i+1]);
      i++;
    }
    if(!string(argv[i]).compare("-M")){
      params.nMutants = atoi(argv[i+1]);
      i++;
    }
    if(!string(argv[i]).compare("-m")){
      params.nPointMutations = atoi(argv[i+1]);
      i++;
    }
    if(!string(argv[i]).compare("-m2")){
      params.nPointMutationsDom2 = atoi(argv[i+1]);
      i++;
    }
    if(!string(argv[i]).compare("-c")){
      params.mutantsListFile = argv[i+1];
      i++;
    }
    if(!string(argv[i]).compare("-ns")){
      params.nDomainSplit = atoi(argv[i+1]);
      i++;
    }
    if(!string(argv[i]).compare("-pi")){
      params.mutatePositionsFile = argv[i+1];
      i++;
    }
    if(!string(argv[i]).compare("-T")){
      params.T = atof(argv[i+1]);
      i++;
    }
    if(!string(argv[i]).compare("-t")){
      params.probThreshold = atof(argv[i+1]);
      i++;
    }
  }  
  this->setOutPrefix(params.outPrefix);
  return params;
}

PottsModel IOHandler::loadModel(string prmFile){
  PottsModel model;

  ifstream prm(prmFile.c_str(),ios::in);
  if(prm){
    // Parse header
    string line;    
    getline(prm,line);
    vector<string> splitLine = splitString(line,' ');
    model.N = atoi(splitLine[1].c_str());
    model.q = atoi(splitLine[2].c_str());
    uInt nParams = model.q*model.N + model.q*model.q*model.N*(model.N-1)/2;

    // Allocate and load parameters
    model.P = new float[nParams];
    for(uInt i=0;i<nParams;i++){
      getline(prm,line);
      splitLine = splitString(line,' ');
      if(i<model.q*model.N)
	model.P[h_at(atoi(splitLine[0].c_str())-1,atoi(splitLine[1].c_str())-1,model.q)] 
	  = atof(splitLine[2].c_str());
      else
	model.P[J_at(atoi(splitLine[0].c_str())-1,atoi(splitLine[1].c_str())-1,
 	       atoi(splitLine[2].c_str())-1,atoi(splitLine[3].c_str())-1,model.N,model.q)]
	  = atof(splitLine[4].c_str());      
    }
  }
  else{
    cout<<"Unable to read prm file"<<endl;
    abort();
  }
  prm.close();  
  return model;
}

Samples IOHandler::loadSamples(string sampleFile){
  ifstream in(sampleFile.c_str());
  string tmp;
  getline(in,tmp);
  vector<string> split = splitString(tmp,' ');
  vector<uInt> raw(0);
  for(uInt i=0;i<split.size();i++)
    raw.push_back(atoi(split[i].c_str()));
  while(!in.eof()){
    in >> tmp;
    raw.push_back(atoi(tmp.c_str()));
  }
  Samples samples;
  samples.N = split.size();
  samples.q = *max_element(raw.begin(),raw.end())+1;
  samples.Nsamples = raw.size()/samples.N;
  samples.data = new uInt[samples.Nsamples*samples.N];
  std::memcpy(samples.data,&raw[0],samples.N*samples.Nsamples*sizeof(uInt));
  
  return samples;
}

void IOHandler::writeConfigurations(uInt* configurations, uInt Nsamples, uInt N){
  ofstream out((outPrefix+"_confs.dat").c_str(),ios::out);
  for(uInt s=0;s<Nsamples;s++){
    for(uInt i=0;i<N;i++)
      out<<configurations[s*N+i]<<" ";
    out<<endl;
  }
  out.close();
}

void IOHandler::writeEnergies(float* Es, uInt Nsamples,uInt nFields){
  ofstream out((outPrefix+"_Es.dat").c_str(),ios::out);
  for(uInt i=0;i<Nsamples;i++){
    for(uInt nF=0;nF<nFields;nF++){
      if(nF<(nFields-1))
	out<<Es[nFields*i+nF]<<"\t";
      else
	out<<Es[nFields*i+nF];
    }
    out<<endl;
  }
  out.close();
}

void IOHandler::writePrm(float* P,uInt N, uInt q){
  ofstream out((outPrefix+".prm").c_str(),ios::out);
  out<<"Header   "<<N<<"   "<<q<<endl;
  for(uInt i=0;i<N;i++)
    for(uInt A=0;A<q;A++)
      out<<"    "<<i+1<<"    "<<A+1<<"  "<<P[h_at(i,A,q)]<<endl;
  for(uInt i=0;i<N-1;i++)
    for(uInt j=i+1;j<N;j++)
      for(uInt A=0;A<q;A++)
	for(uInt B=0;B<q;B++)
	  out<<"    "<<i+1<<"   "<<j+1<<"   "<<A+1<<"   "<<B+1<<"   "<<P[J_at(i,j,A,B,N,q)]<<endl;
  out.close();
}

void IOHandler::writeFrequencies(float* F,uInt N, uInt q){
  ofstream out((outPrefix+"_fi.dat").c_str(),ios::out);
  for(uInt i=0;i<N;i++)
    for(uInt A=0;A<q;A++)
      out<<F[h_at(i,A,q)]<<endl;
  out.close();
  out.open((outPrefix+"_fij.dat").c_str(),ios::out);
  for(uInt i=0;i<N-1;i++)
    for(uInt j=i+1;j<N;j++)
      for(uInt A=0;A<q;A++)
	for(uInt B=0;B<q;B++)
	  out<<F[J_at(i,j,A,B,N,q)]<<endl;
  out.close();
}
