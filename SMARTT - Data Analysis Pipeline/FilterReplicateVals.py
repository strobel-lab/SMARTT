# -*- coding: utf-8 -*-
"""
Python script for comparing the paramter values output for two replicates and 
ouputing an excel file with three tabs:
    - Tab 1 contains the unfiltered data for both replicates
    - Tab 2 filters the data if the error value for either K_1/2 or amplitude 
      are larger than the associated parameter value.
    - Tab 3 provides the difference map (shows the filtered values)
 
Running this script:
 - CSV files for the two replicates
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

#import math
#import operator
import os
import pandas as pd
import numpy as np


# File locations
rep1Folder = os.path.join("demos", "Bsu_data")
rep2Folder = os.path.join("demos", "Bsu_data")
outFolder = "output"

rep1File = os.path.join(rep1Folder, "BsuR1-WT_QST20_Fits_SinglesOnly.csv")
rep2File = os.path.join(rep2Folder, "BsuR2-WT_QST20_Fits_SinglesOnly.csv")
outFile = os.path.join(outFolder, "Bsu-WT_RepValsAnalysis.xlsx")

#rep1File = os.path.join(rep1Folder, "BsuR1-U87G_QST20_Fits_SinglesOnly.csv")
#rep2File = os.path.join(rep2Folder, "BsuR2-U87G_QST20_Fits_SinglesOnly.csv")
#outFile = os.path.join(outFolder, "Bsu-U87G_RepValsAnalysis.xlsx")

#rep1File = os.path.join(rep1Folder, "BsuR1-U185G_QST20_Fits_SinglesOnly.csv")
#rep2File = os.path.join(rep2Folder, "BsuR2-U185G_QST20_Fits_SinglesOnly.csv")
#outFile = os.path.join(outFolder, "Bsu-U185G_RepValsAnalysis.xlsx")

# Read input files 
R1_df = pd.read_csv(rep1File, skipinitialspace=True, header=0, index_col = 0,
                    keep_default_na=False, na_values=[""])

R2_df = pd.read_csv(rep2File, skipinitialspace=True, header=0, index_col = 0,
                      keep_default_na=False, na_values=[""])

#Initialize Values
R = 0.001987 # kcal / (K mol)
T = 310.15   # K

# Calculate the apparent free energy based on K
R1_df["dG_app"] = R*T*np.log(R1_df['K'])
R1_df["-"] = np.nan     #empty column to make visualization of the Rep1 and Rep2 data better in the output file
R2_df["dG_app"] = R*T*np.log(R2_df['K'])

# Create the unfiltered dataframe by merging the two datasets (type = outer)
comb_df1 = pd.merge(R1_df.add_prefix('R1_'), R2_df.add_prefix('R2_'), how='outer', left_index=True, right_index=True)

# Copy dataframes for making the filtered values and difference map sheets
comb_df2 = comb_df1.copy()
comb_df3 = comb_df1.copy()

## Filter values that have poor fits
#R1_df.loc[(R1_df['K'] < R1_df['K_SD']) | (R1_df['Amp'] < R1_df['Amp_SD'])] = np.nan # Note that '|' must be used here, not 'or'
#R2_df.loc[(R2_df['K'] < R2_df['K_SD']) | (R2_df['Amp'] < R2_df['Amp_SD'])] = np.nan # Note that '|' must be used here, not 'or'

# Filter values if either fit was bad
comb_df2.loc[(comb_df1['R1_K'] < comb_df1['R1_K_SD']) | (comb_df1['R1_Amp'] < comb_df1['R1_Amp_SD'])
| (comb_df1['R2_K'] < comb_df1['R2_K_SD']) | (comb_df1['R2_Amp'] < comb_df1['R2_Amp_SD'])] = np.nan

# Filter values that have good fits (for creating a difference map)
comb_df3.loc[(comb_df1['R1_K'] > comb_df1['R1_K_SD']) & (comb_df1['R1_Amp'] > comb_df1['R1_Amp_SD'])
& (comb_df1['R2_K'] > comb_df1['R2_K_SD']) & (comb_df1['R2_Amp'] > comb_df1['R2_Amp_SD'])] = np.nan

# Create a Pandas Excel writer using XlsxWriter as the engine.
writer = pd.ExcelWriter(outFile, engine='xlsxwriter')

# Write each dataframe to a new sheet
comb_df1.to_excel(writer, sheet_name='Unfiltered')
comb_df2.to_excel(writer, sheet_name='Filtered')
comb_df3.to_excel(writer, sheet_name='Difference Map') #Useful for displaying which variants had 1 or more fits filtered

# Write outfile and close the excel writer
writer.save()