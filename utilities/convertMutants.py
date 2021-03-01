#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" C3D: Computational protein Design by Duplication and Divergence (2021)
    Author: Duccio Malinverni, St.Jude Children's Research Hospital, Memphis, TN, USA 
    This file is covered by GPL-3.0 license (see LICENSE file in the root of this project.    
"""

import sys
import numpy as np
import copy

if  __name__=='__main__':

    mutantsFile=sys.argv[1]
    nativeFile=sys.argv[2]
    outPrefix=sys.argv[3]

    aaDict={0:'-', 1:'A', 2:'C', 3:'D', 4:'E', 5:'F', 6:'G', 7:'H', 8:'I', 9:'K', 10:'L', 11:'M',
            12:'N', 13:'P', 14:'Q', 15:'R', 16:'S', 17:'T', 18:'V', 19:'W', 20:'Y'}   

    native=np.loadtxt(nativeFile)
    fasta=open(outPrefix+'_mutants.fasta','w')
    mutants=open(outPrefix+'_mutantsAA.dat','w')

    with open(mutantsFile,'r') as f:
        for il,line in enumerate(f):
            if line[0]=="#":
                mutants.write(line)
            else:
                ls=line.strip().split('\t')
                mutations=ls[:-5]
                fasta.write('>'+str(il)+'\n')
                seq=copy.copy(native)
                for m in mutations:
                    ms=m.split('_')
                    mutants.write(aaDict[int(ms[0])]+str(1+int(ms[1]))+aaDict[int(ms[2])]+'\t')
                    seq[int(ms[1])]=int(ms[2])
                mutants.write('\t'.join(ls[-5:])+'\n')
                fasta.write(''.join(aaDict[s] for s in seq)+'\n')
                
    fasta.close()
    mutants.close()
