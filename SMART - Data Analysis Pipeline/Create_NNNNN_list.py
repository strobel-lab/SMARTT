
# Python script for analyzing the output from Create_NNNNN_List.py to create a csv 
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
#import re
#import ast


NUCS = ["A", "T", "G", "C"]
COMP = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C"
}

def reverse_complement(seq):
    # Complement sequence
    comp = "".join([COMP.get(c,"N") for c in seq])
    RC = comp[::-1]
    return RC

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~Variables that can be edited~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# If debugging, set to True and file inputs will be ignored
IgnoreFileInputs = True

# These must be set if IgnoreFileInputs is set to True
# Files and input data
if IgnoreFileInputs == True:
    demoDirectory = "demos"
    fastqList = [os.path.join(demoDirectory, "Cte_100mM_Gly_demo.fastq")]
    outDirectory = "output"
    outputFileName = os.path.join(outDirectory, "Cte_100mM_Gly_demo_NNNNN_list.txt")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

# Check file paths inputs
if len(sys.argv) != 3:
    InputLength = False
else:
    InputLength = True

if InputLength == False and IgnoreFileInputs == False:
    sys.exit("Usage: python Create_NNNNN_list.py <fastq file(s)> <output filename>")
else:          
    # Files and input data
    if IgnoreFileInputs == False:
        fastqList = sys.argv[1].split(',')
        outputFileName = sys.argv[2]

# Determine behavior if output file already exists/open the outfile
if os.path.isfile(outputFileName):
    sys.stderr.write("Output file already exists; overwriting file!\n\n")
    fileOut = open(outputFileName, "w")
else:
    fileOut = open(outputFileName, "w")

# Initialize variables
SeqLine = False
QNAME = ''
NNNNN_identity = ''
line2write = ''



for fastqPath in fastqList:
    
    # Read and analyze each specific file in the file list
    with open(fastqPath, 'r') as f:
        
        # Determine QName and NNNNN values
        for line in f:    
        # Skip headers
            if line[0] == "@":
                parsedLine = line.split()
                QNAME = str(parsedLine[0][1:])
                SeqLine = True
                continue
            
            if SeqLine == True:
                NNNNN_identity = str(reverse_complement(line[:5]))
                output = '\t'.join([QNAME, NNNNN_identity]) + '\n'
                
                # Write NNNNN info to file
                fileOut.write(output)
                
                # Reinitialize variables
                output = ''
                QNAME = ''
                NNNNN_identity = ''
                SeqLine = False
            
    
            
fileOut.close()
    