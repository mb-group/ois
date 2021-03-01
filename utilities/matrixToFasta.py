#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" C3D: Computational protein Design by Duplication and Divergence (2021)
    Author: Duccio Malinverni, St.Jude Children's Research Hospital, Memphis, TN, USA 
    This file is covered by GPL-3.0 license (see LICENSE file in the root of this project.    
"""

import sys
import Bio.SeqIO
import Bio.SeqRecord
import numpy as np

if  __name__=='__main__':

    matrixFile = sys.argv[1]
    fastaFile = sys.argv[2]

    aaDict={0:'-', 1:'A', 2:'C', 3:'D', 4:'E', 5:'F', 6:'G', 7:'H', 8:'I', 9:'K', 10:'L', 11:'M',
            12:'N', 13:'P', 14:'Q', 15:'R', 16:'S', 17:'T', 18:'V', 19:'W', 20:'Y'}   

    msa=[]
    with open(matrixFile,'r') as f:
        for i,line in enumerate(f):
            seq =  ''.join([aaDict[int(aa)] for aa in line.split()])
            msa.append(Bio.SeqRecord.SeqRecord(Bio.Seq.Seq(seq),id=str(i),name=str(i),description=str(i)))

    Bio.SeqIO.write(msa,fastaFile,"fasta")

    
