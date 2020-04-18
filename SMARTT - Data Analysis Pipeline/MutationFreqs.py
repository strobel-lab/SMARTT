# -*- coding: utf-8 -*-
"""
This python script takes a Percent Full Length csv file for a single replicate
and prints the number of AT & GC base pair positions as well as the observed 
mutation frequency. It requires that only the single mutations are contained 
within the file.

Running this script:
 - Requires the Percent Full Length csv file for a single replicate (single 
   mutations only)
 - Requires Python 3 (has been tested with version 3.6)

Notes:
 - Public release 1.0
 - Copyright 2020 Chad Torgerson

##################################################################################
 GPL statement:

 This file is part of SMARTT-DAP (SMARTT - Data Analysis Pipeline).
 
 SMARTT-DAP is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 SMARTT-DAP is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
##################################################################################
"""

import re
import os
import numpy as np
import pandas as pd

try:
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~Variables to edit~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    # File locations
    folder = os.path.join("demos", "Bsu_data")
    variantFile ="BsuR1-WT_QST20_PFL_SinglesOnly.csv"
#    variantFile ="BsuR2-WT_QST20_PFL_SinglesOnly.csv"
#    variantFile ="BsuR1-U87G_QST20_PFL_SinglesOnly.csv"
#    variantFile ="BsuR2-U87G_QST20_PFL_SinglesOnly.csv"
#    variantFile ="BsuR1-U185G_QST20_PFL_SinglesOnly.csv"
#    variantFile ="BsuR2-U185G_QST20_PFL_SinglesOnly.csv"
    variantPath = os.path.join(folder, variantFile)
    

    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    # Reading variant files
    df=pd.read_csv(variantPath, sep=',',header=0)
       
    # Initializing variables
    mutDict = {}
    nucCount = {}
    mutList = [
            'A-->G, T-->C',
            'A-->T, T-->A',
            'A-->C, T-->G',
            'G-->A, C-->T',
            'G-->T, C-->A',
            'G-->C, C-->G',
            ]
    
    # Fill in first column with the Sample Concentrations
    sampleConcList = list(df.values[0,2::2])
    outputMatrix = np.zeros((len(sampleConcList),len(df)*2-1))
    outputMatrix[:,0] = sampleConcList
    
    count = 0
    for df_row in range(1,len(df)):
        varName = df.values[df_row][0]
        if varName == 'WT':
            continue
#        headers.extend([varName])
        varMatch = re.match(r"([a|u|t|c|g|A|U|T|C|G])[0-9]+([a|u|t|c|g|A|U|T|C|G])", varName)
        mutType = (varMatch.group(1),varMatch.group(2))
        nucType = varMatch.group(1).upper()
        
        sampleReads = df.values[df_row][1::2]
        sampleReads = np.array([float(i) for i in sampleReads])
        total = np.sum(sampleReads)
        
        try:
            mutDict[mutType] += total
        except:
            mutDict[mutType] = total
            
        try:
            nucCount[nucType] += 1
        except:
            nucCount[nucType] = 1
        
    nucCount = {k: v / 3 for k, v in nucCount.items()}
    AT = nucCount['A'] + nucCount['T']
    GC = nucCount['G'] + nucCount['C']
    print("AT base pair positions:", str(int(AT)))
    print("GC base pair positions:", str(int(GC)))
    
    mutDict2 = {
            'A-->G, T-->C' : (mutDict[('A', 'G')] + mutDict[('T', 'C')])/AT,
            'A-->T, T-->A' : (mutDict[('A', 'T')] + mutDict[('T', 'A')])/AT,
            'A-->C, T-->G' : (mutDict[('A', 'C')] + mutDict[('T', 'G')])/AT,
            'G-->A, C-->T' : (mutDict[('G', 'A')] + mutDict[('C', 'T')])/GC,
            'G-->T, C-->A' : (mutDict[('G', 'T')] + mutDict[('C', 'A')])/GC,
            'G-->C, C-->G' : (mutDict[('G', 'C')] + mutDict[('C', 'G')])/GC,
            }
    totalCounts = sum(mutDict2.values())
    mutationFreqs = {k: 100 * v / totalCounts for k, v in mutDict2.items()}
    
    for item in mutList:
        print(item, round(mutationFreqs[item],2))
    
    
except:
    raise