#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" C3D: Computational protein Design by Duplication and Divergence (2021)
    Author: Duccio Malinverni, St.Jude Children's Research Hospital, Memphis, TN, USA 
    This file is covered by GPL-3.0 license (see LICENSE file in the root of this project.    
"""

import sys
import Bio.SeqIO
import numpy as np

if  __name__=='__main__':

    fastaFile = sys.argv[1]
    matrixFile = sys.argv[2]
    msa=Bio.SeqIO.parse(fastaFile,'fasta')
        
    seqs=[list(str(s.seq)) for s in msa]     
    aaDict={'-':0, 'X':0,'Z':0,'B':0, 'A':1, 'C':2, 'D':3, 'E':4, 'F':5, 'G':6, 'H':7, 'I':8, 'K':9, 'L':10, 'M':11,
            'N':12, 'P':13, 'Q':14, 'R':15, 'S':16, 'T':17, 'V':18, 'W':19, 'Y':20}   

    matrix= np.asarray([[aaDict[aa] for aa in seq] for seq in seqs])
    np.savetxt(matrixFile,matrix,fmt='%d',delimiter=' ')
    
