# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 11:05:11 2020

@author: justi
"""

def BasicSpars(WingShape,SparLoc):
    import numpy as np
    nSpan  = len(WingShape)
    nSpars = len(SparLoc)
    # WingStruct = np.zeros((nSpan,nSpars+1))
    WingStruct = [''] * nSpan
    
    for i in range(nSpan):
        if np.size(SparLoc) == nSpars:  # for constant spar locations
            WingStruct[i] = SparLoc
        else:  # for variable spar locations
            WingStruct[i] = SparLoc[i]
    return WingStruct
        
        
