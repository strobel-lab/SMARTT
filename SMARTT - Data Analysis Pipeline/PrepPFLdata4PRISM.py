
# Python script for preparing the output from DeterminePFL.py for analysis in PRISM
# with the end goal of fitting the transcription termination data  
# (Step 3 in the pipeline)
#
# - Requires Python 3 (has been tested with version 3.5)
# - Public release 1.0
# - Copyright 2018 Chad Torgerson

###################################################################################
# GPL statement:
#
# This file is part of SMARTT-DAP (SMARTT - Data Analysis Pipeline).
# 
# SMARTT-DAP is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SMARTT-DAP is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
###################################################################################


import os
import numpy as np
import pandas as pd

try:
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~Variables to edit~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    # File locations
    folder = 'demos'
    variantFile ="Cte_PFL-demo.csv"
    variantPath = os.path.join(folder, variantFile)
    
    outputDirectory = 'output'
    sampleName = "Cte_PFL-demo"
    outputFileName = os.path.join(outputDirectory, sampleName + "_prepped4PRISM.xlsx")
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    # Reading variant files
    df=pd.read_csv(variantPath, sep=',',header=0)
       
    # Initializing header
    headers = ['Conc']
    
    # Fill in first column with the Sample Concentrations
    sampleConcList = list(df.values[0,2::2])
    outputMatrix = np.zeros((len(sampleConcList),len(df)*2-1))
    outputMatrix[:,0] = sampleConcList
    
    
    for df_row in range(1,len(df)):
        varName = df.values[df_row][0]
        if varName != 'WT':
            varName = varName.replace('T','U')
        headers.extend([varName,varName+'_stdev'])
        
        sampleReads = df.values[df_row][1::2]
        samplePFLs = df.values[df_row][2::2]
        sampleReads = np.array([float(i) for i in sampleReads])
        samplePFLs = np.array([float(i)/100 for i in samplePFLs])
        isapprop_1 = np.nanmin(samplePFLs*sampleReads)
        isapprop_2 = np.nanmin((1-samplePFLs)*sampleReads)
        
        # Checks to see if n*p < 5 or n(1-p) < 5      <--- Criteria to decide if sample size is too small
        # Bayesian inference is used for values that don't have a sufficient number of reads
        if max(sampleReads) == 0:
            continue
        elif isapprop_1 < 5 or isapprop_2 < 5:
            successes = np.around(samplePFLs*sampleReads)
            samplePFLs = (successes + 0.5)/(sampleReads + 1)
            PFL_stdev = np.sqrt(((successes + 0.5)*(sampleReads-successes + 0.5))/(((sampleReads + 1)**2)*(sampleReads + 2)))
        else:
            PFL_stdev = np.sqrt(samplePFLs*(1-samplePFLs)/sampleReads); # Standard Deviation of the probability. This is equal to -> sqrt(p*(1-p)/n)
        
        # Fill in the PFLs and stdev
        outputMatrix[:,2*df_row-1] = samplePFLs*100
        outputMatrix[:,2*df_row] = PFL_stdev*100
    
    output_df = pd.DataFrame(outputMatrix) # outputMatrix is a numpy 2d array
    writer = pd.ExcelWriter(outputFileName, engine='xlsxwriter')
    output_df.to_excel(writer, header=headers,index=False) # headers is a list of string corresponding to the title of each column of output_df
    writer.save()
    
except:
    raise