# -*- coding: utf-8 -*-
"""
Created on Wed Aug 22 15:31:58 2018

@author: Chad
"""


# Python script for analyzing the output from Add_NNNNN_Vals2VariantList.py to create a csv 
# file with the number of reads and the percent of full length RNA produced by a 
# library of mutant RNA (Step 2 in the Pipeline)
#
# - Requires Python 3 (has been tested with version 3.5)
# - Public release 1.0
# - Copyright 2018 Chad Torgerson

###################################################################################
# GPL statement:
#
# This file is part of SMART-DAP (SMART - Data Analysis Pipeline).
# 
# SMART-DAP is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SMART-DAP is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
###################################################################################

import sys
import os

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~Variables that can be edited~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# If debugging, set to True and file inputs will be ignored
IgnoreFileInputs = False

# These must be set if IgnoreFileInputs is set to True
# Files and input data
if IgnoreFileInputs == True:
    demoDirectory = "demos"
    variant_inFile = os.path.join(demoDirectory, "Cte_100mM_Gly-demo_variants.txt")
    NNNNN_inFile = os.path.join(demoDirectory, "Cte_100mM_Gly_NNNNN_list_test.txt")
    outDirectory = "output"
    outputFileName = os.path.join(outDirectory, "Cte_100mM_Gly-demo_variants_NNNNNs-added.txt")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 


# Check file paths inputs
if len(sys.argv) != 4:
    InputLength = False
else:
    InputLength = True

if InputLength == False and IgnoreFileInputs == False:
    sys.exit("Usage: python Add_NNNNN_Vals2VariantList.py <variant list txt file> <NNNNN_list txt file> <output filename>")
else:          
    # Files and input data
    if IgnoreFileInputs == False:
        variant_inFile = sys.argv[1]
        NNNNN_inFile = sys.argv[2]
        outputFileName = sys.argv[3]
        print(variant_inFile)
        print(NNNNN_inFile)
        print(outputFileName)

# Save sample name; remove path and filename extension
sampleName = os.path.splitext(os.path.split(variant_inFile)[1])[0]

# Initialize variables
NNNNN_dict = {}
count = 0
XXXXX_count = 0
QNAME = 'Qval'
NNNNN_identity = 'XXXXX'
newLine = ''
header = "\t".join([
            "variant",
            "(start, stop)",
            "QNAME",
            "NNNNN_value"+ "\n"
            ])

    
with open(NNNNN_inFile, 'r') as f1:
    for line in f1:
       (key, val) = line.split()
       NNNNN_dict[key] = val

# Determine behavior if output file already exists/open the outfile
if os.path.isfile(outputFileName):
    sys.stderr.write("Output file already exists; overwriting file!\n\n")
    fileOut = open(outputFileName, "w")
else:
    fileOut = open(outputFileName, "w")

with open(variant_inFile, 'r') as f2:
    for line in f2:
        if count == 0:
            count += 1
            # Write header to output file
            fileOut.write(header)
            continue
        
        # Parse Line
        parsedLine = line.split('\t')
        
        # Determine QNAME (for lookup in NNNNN_dict)
        QNAME = parsedLine[2][:-1]
        
        # Look up NNNNN value; if not found, label as 'XXXXX'
        try:
            NNNNN_identity = NNNNN_dict[QNAME]
        except:
            NNNNN_identity = 'XXXXX'
            XXXXX_count += 1
        
        # Write output (with NNNNN value appended)
        newLineList = parsedLine[:2]
        newLineList.extend((QNAME, NNNNN_identity))
        newLine = '\t'.join(newLineList) + '\n'
        fileOut.write(newLine)

fileOut.close()
print('File:', sampleName)
print('XXXXX Count:', XXXXX_count, '\n')