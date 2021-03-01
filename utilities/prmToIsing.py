#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" C3D: Computational protein Design by Duplication and Divergence (2021)
    Author: Duccio Malinverni, St.Jude Children's Research Hospital, Memphis, TN, USA 
    This file is covered by GPL-3.0 license (see LICENSE file in the root of this project.    
"""

import sys
import numpy as np
import pandas as pd

if  __name__=='__main__':

    inPrm=sys.argv[1]
    outPrm=sys.argv[2]
    
    # Load prm file
    with open(inPrm,'r') as f:
        line = f.readline()
        splt = line.split()
        N = int(splt[1])
        q = int(splt[2])
        
        with open(inPrm,'r') as f:
            h = pd.read_csv(f,usecols=[2],delimiter='\s+',nrows=N*q,engine='python').values
        with open(inPrm,'r') as f:
            J = pd.read_csv(f,usecols=[2],delimiter='\s+',nrows=N*(N-1)*q*q,
                            skiprows=N*q,engine='python').values
        h=h.ravel()
        J=J.ravel()

    h0 = np.zeros(h.shape)
    J0 = np.zeros(J.shape)

    # Shift couplings to zero-sum gauge
    c=0
    for i in range(N):
        h0[i*q:(i+1)*q] += h[i*q:(i+1)*q]
        for j in range(i+1,N):
            jij = np.reshape(J[c:c+q*q],(q,q))
            j0 = jij - jij.mean(axis=0,keepdims=True) - jij.mean(axis=1,keepdims=True) + jij.mean()
            J0[c:c+q*q]=np.reshape(j0,(q*q))
            c+=q*q
            
            h0[i*q:(i+1)*q] += jij.mean(axis=1)
            h0[j*q:(j+1)*q] += jij.mean(axis=0)
            
    # Shift fields to zero-sum gauge
    for i in range(N):
        h0[i*q:(i+1)*q]-=h0[i*q:(i+1)*q].mean()

    # Save shifted parameters
    with open(outPrm,'w') as f:
        f.write('protein   '+str(N)+'   '+str(q)+'\n')
        for i in range(N):
            for A in range(q):
                f.write('    '+str(i+1)+'    '+str(A+1)+'    '+str(h0[q*i+A])+'\n')
        c=0
        for i in range(N):
            for j in range(i+1,N):
                for A in range(q):
                    for B in range(q):
                        f.write('    '+str(i+1)+'    '+str(j+1)+'    '+str(A+1)+'    '+str(B+1)
                                +'    '+str(J0[c])+'\n')
                        c+=1

