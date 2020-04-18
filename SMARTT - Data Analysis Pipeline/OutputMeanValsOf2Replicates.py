# -*- coding: utf-8 -*-
"""
This python script for compares the paramter values output for two replicates and 
ouputs a CSV file containing the mean parameter value and error values across the
two replicates. For samples containing error values larger than the parameter value
for amplitude or K, only Ymin is retained (the other parameters are filtered).
For samples with a Hill coefficient, the Hill coefficent may have been filtered 
(even if the other parameters were not) if the error term was larger than the 
Hill value or if the error term was greater than 1.
 
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

rep1File = os.path.join(rep1Folder, "BsuR1-WT_Fits_Hill-SinglesOnly.csv")
rep2File = os.path.join(rep2Folder, "BsuR2-WT_Fits_Hill-SinglesOnly.csv")
outFile = os.path.join(outFolder, "Bsu-WT_MeanFitValues-withHill-R.csv")
hasHill = True

#rep1File = os.path.join(rep1Folder, "BsuR1-U185G_QST20_Fits_SinglesOnly.csv")
#rep2File = os.path.join(rep2Folder, "BsuR2-U185G_QST20_Fits_SinglesOnly.csv")
#outFile = os.path.join(outFolder, "Bsu-U185G_MeanFitValues-R.csv")
#hasHill = False

#rep1File = os.path.join(rep1Folder, "BsuR1-U87G_QST20_Fits_SinglesOnly.csv")
#rep2File = os.path.join(rep2Folder, "BsuR2-U87G_QST20_Fits_SinglesOnly.csv")
#outFile = os.path.join(outFolder, "Bsu-U87G_MeanFitValues.csv")
#hasHill = False

# Read input files 
R1_df = pd.read_csv(rep1File, skipinitialspace=True, header=0, index_col = 0,
                    keep_default_na=False, na_values=[""])

R2_df = pd.read_csv(rep2File, skipinitialspace=True, header=0, index_col = 0,
                      keep_default_na=False, na_values=[""])

R1_df = R1_df.apply(pd.to_numeric, errors='coerce')
R2_df = R2_df.apply(pd.to_numeric, errors='coerce')

## Copy column names for later
#headers = R1_df.columns

#R1_df.loc[(R1_df['K'] < R1_df['K_SD']) | (R1_df['Amp'] < R1_df['Amp_SD']), ~R1_df.columns.isin(['Ymin', 'Ymin_SD'])] = np.nan # Note that '|' must be used here, not 'or'
R1_df['Reps'] = 1
R2_df['Reps'] = 1
R1_df.loc[(R1_df['K'] < R1_df['K_SD']) | (R1_df['Amp'] < R1_df['Amp_SD']), ['Reps', 'K', 'K_SD', 'Ymax', 'Ymax_SD', 'Amp', 'Amp_SD']] = np.nan # Note that '|' must be used here, not 'or'
R2_df.loc[(R2_df['K'] < R2_df['K_SD']) | (R2_df['Amp'] < R2_df['Amp_SD']), ['Reps', 'K', 'K_SD', 'Ymax', 'Ymax_SD', 'Amp', 'Amp_SD']] = np.nan # Note that '|' must be used here, not 'or'

if hasHill == True:
    R1_df.loc[(R1_df['K'] < R1_df['K_SD']) | (R1_df['Amp'] < R1_df['Amp_SD']), ['Reps', 'K', 'K_SD', 'Ymax', 'Ymax_SD', 'Amp', 'Amp_SD', 'Hill', 'Hill_SD']] = np.nan # Note that '|' must be used here, not 'or'
    R2_df.loc[(R2_df['K'] < R2_df['K_SD']) | (R2_df['Amp'] < R2_df['Amp_SD']), ['Reps', 'K', 'K_SD', 'Ymax', 'Ymax_SD', 'Amp', 'Amp_SD', 'Hill', 'Hill_SD']] = np.nan # Note that '|' must be used here, not 'or'

    R1_df['Hill_Reps'] = 1
    R2_df['Hill_Reps'] = 1
    R1_df.loc[(R1_df['Hill'] < R1_df['Hill_SD']) | (R1_df['Hill_SD'] > 1), ['Hill_Reps', 'Hill', 'Hill_SD']] = np.nan 
    R2_df.loc[(R2_df['Hill'] < R2_df['Hill_SD']) | (R2_df['Hill_SD'] > 1), ['Hill_Reps', 'Hill', 'Hill_SD']] = np.nan
else:
    R1_df.loc[(R1_df['K'] < R1_df['K_SD']) | (R1_df['Amp'] < R1_df['Amp_SD']), ['Reps', 'K', 'K_SD', 'Ymax', 'Ymax_SD', 'Amp', 'Amp_SD']] = np.nan # Note that '|' must be used here, not 'or'
    R2_df.loc[(R2_df['K'] < R2_df['K_SD']) | (R2_df['Amp'] < R2_df['Amp_SD']), ['Reps', 'K', 'K_SD', 'Ymax', 'Ymax_SD', 'Amp', 'Amp_SD']] = np.nan # Note that '|' must be used here, not 'or'


    
# Add a prefix to differentiate the two datasets that are about to be merged
R1_df = R1_df.add_prefix('R1_')
R2_df = R2_df.add_prefix('R2_')

#R1_df.columns = headers

# Merge the two datasets (type = outer)
comb_df = pd.merge(R1_df, R2_df, how='outer', left_index=True, right_index=True)

# Initialize final dataframe
if hasHill == False:
    cols = ['Reps', 'K', 'K_SD', 'Ymin', 'Ymin_SD', 'Ymax', 'Ymax_SD', 'Amp', 'Amp_SD']
else:
    cols = ['Reps', 'K', 'K_SD', 'Ymin', 'Ymin_SD', 'Ymax', 'Ymax_SD', 'Amp', 'Amp_SD', 'Hill_Reps', 'Hill', 'Hill_SD']
inds = comb_df.index
final_df = pd.DataFrame(index = inds, columns = cols)

# Fill final dataframe
final_df['Reps'] = comb_df[['R1_Reps', 'R2_Reps']].sum(axis=1)
final_df['K'] = comb_df[['R1_K', 'R2_K']].mean(axis=1)
final_df['K_SD'] = comb_df[['R1_K', 'R2_K']].std(axis=1, ddof=1)
final_df['Ymin'] = comb_df[['R1_Ymin', 'R2_Ymin']].mean(axis=1)
final_df['Ymin_SD'] = comb_df[['R1_Ymin', 'R2_Ymin']].std(axis=1, ddof=1)
final_df['Ymax'] = comb_df[['R1_Ymax', 'R2_Ymax']].mean(axis=1)
final_df['Ymax_SD'] = comb_df[['R1_Ymax', 'R2_Ymax']].std(axis=1, ddof=1)
final_df['Amp'] = comb_df[['R1_Amp', 'R2_Amp']].mean(axis=1)
final_df['Amp_SD'] = comb_df[['R1_Amp', 'R2_Amp']].std(axis=1, ddof=1)

if hasHill == True:
    final_df['Hill_Reps'] = comb_df[['R1_Hill_Reps', 'R2_Hill_Reps']].sum(axis=1)
    final_df['Hill'] = comb_df[['R1_Hill', 'R2_Hill']].mean(axis=1)
    final_df['Hill_SD'] = comb_df[['R1_Hill', 'R2_Hill']].std(axis=1, ddof=1)

# If n = 1 
final_df.loc[(final_df.Reps == 1), 'K_SD'] = comb_df[['R1_K_SD', 'R2_K_SD']].mean(axis=1)
final_df.loc[(final_df.Reps == 1), 'Ymin_SD'] = comb_df[['R1_Ymin_SD', 'R2_Ymin_SD']].mean(axis=1)
final_df.loc[(final_df.Reps == 1), 'Ymax_SD'] = comb_df[['R1_Ymax_SD', 'R2_Ymax_SD']].mean(axis=1)
final_df.loc[(final_df.Reps == 1), 'Amp_SD'] = comb_df[['R1_Amp_SD', 'R2_Amp_SD']].mean(axis=1)
if hasHill == True:
    final_df.loc[(final_df.Hill_Reps == 1), 'Hill_SD'] = comb_df[['R1_Hill_SD', 'R2_Hill_SD']].mean(axis=1)


final_df.to_csv(outFile)