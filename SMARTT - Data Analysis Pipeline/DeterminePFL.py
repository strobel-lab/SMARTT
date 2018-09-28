
# Python script for analyzing the output from CreateVariantsList.py to create a csv 
# file with the number of reads and the percent of full length RNA produced by a 
# library of mutant RNA (Step 2 in the Pipeline)
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



import sys
import os
import re
import ast


try:
    
    def conc2MolValue(concValue, units):
        # Convert to molar concentration and return numerical value (without units)
        if units == 'M':
            MolarValue = concValue
        elif units == 'mM':
            MolarValue = concValue*(10**-3)
        elif units == 'uM':
            MolarValue = concValue*(10**-6)
        elif units == 'nM':
            MolarValue = concValue*(10**-9)
        elif units == 'pM':
            MolarValue = concValue*(10**-12)
        elif units == 'fM':
            MolarValue = concValue*(10**-15)
        return MolarValue    
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~Variables to edit~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    # File locations
    folder = 'demos'
    outputDirectory = 'output'
    variantFileList = (['Cte_100mM_Gly-demo_variants.txt'])
    #variantFileList = (['Cte_100mM_Gly_variants_withNNNNNs_demo.txt'])
    
    #sampleName = "Cte-demo" # This is used to determine the filenames used for keeping track of the sites of termination and the overall percentage of full-length RNA
    sampleName = "Cte_TTS25_QST20"
    outputFileName = os.path.join(outputDirectory, sampleName + "_PFL.csv")
    
    # If file names don't contain the sample concentration, 
    # manually add all of them here and change 'filesContainConc' to False
    # ^This option was never finished and currently doesn't work
    sampleConcList = []
    filesContainConc = True
 
    # Positions for determining if a read should be counted
    begin = 25  # 1st nt position must be always <= this number (this number = first mutation we care about)
    trunc = 190 # last nt must be >= this number and < FL	(190 - 25 = 165)
    FL = 199    # last nt must be >= this number			(199 - 25 = 174)
    offset = 25
    length = 248 # Last possible nt of riboswitch
    
    # Maximum number of mutations to consider
    mutationsAllowed = 1
    
    # If set to true, files will be output reporting the NNNNN values by variant type
    ReportNNNNNs = False
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    # Define pattern for what a concentration looks like
    concentrationRegRE = re.compile(r"([0-9]+)[^a-zA-Z]??([m|u|n|p|f]??M)")
    
    #%%    
    # Initializing variables
    variantSet = set() # Defining an empty set that will keep track of all variant names    
    ReadsPFL_Dict = {} # Prepare Data Dict
    headers = ['variant'] # Assign first header variable as 'variant'
    conc_list = ['Conc.'] # The second line of the file will keep track of the sample concentration values
    LastntHeaders = ['variant'] # Assign first header variable as 'variant' (for sep. file that keeps track of the last positions)
    NNNNN_valsList = []
    NNNNN_indexDict = {}
    Ns_counter = 0
    NUCS = ["A", "C", "G", "T"]
    
    positionNums = list(range(-offset,length+1-offset))
    positionNums.remove(0)
    LastntHeaders.extend(positionNums)
    
    
    if ReportNNNNNs == True:
        for n1 in NUCS:
            for n2 in NUCS:
                for n3 in NUCS:
                    for n4 in NUCS:
                        for n5 in NUCS:
                            NNNNN = n1 + n2 + n3 + n4 + n5
                            NNNNN_valsList.append(NNNNN)
                            NNNNN_indexDict[NNNNN] = Ns_counter
                            Ns_counter += 1
    NNNNN_valsList.append('XXXXX')
    NNNNN_indexDict['XXXXX'] = Ns_counter
    
    # Assign first header variable as 'variant' (for sep. files that keeps track of the NNNNN values)
    NNNNN_Headers = ['variant'] + NNNNN_valsList
    
    
    #%%
    # This is where my loop for iterating through all files starts
    for file in variantFileList:
        sumReadsParsed = 0
        sumReadsAnalyzed = 0
        sumBeganLate = 0
        sumEarlyTrunc = 0
        zeroMutations = 0
        oneMutation = 0
        twoMutations = 0
        threeplusMutations = 0
        LastntDict = {}
        NNNNN_Tnc_Dict = {}
        NNNNN_FL_Dict = {}
        NNNNN_Combined_Dict = {}
        
        
        # Make a list of all the sample concentrations
        # (Searches current file name for concentration and appends it to the list)      
        if filesContainConc == True:
            conc = concentrationRegRE.search(file)            
            try:
                sampleConc = conc.group(1)+conc.group(2)
                if "gly" in file.casefold():
                    sampleConc += "_G"
                elif "ala" in file.casefold():
                    sampleConc += "_A"
                sampleConcList.append(sampleConc)

                try:
                    molarConcValue = conc2MolValue(int(conc.group(1)), conc.group(2))
                except:
                    molarConcValue = 'NA'
            except:
                print("File does not contain concentration in name")
                sampleConcList.append('index' + str(variantFileList.index(file)))
        else:
            sampleConc = sampleConcList[variantFileList.index(file)]
            molarConcValue = 'NA'
        
        # Setup Truncation Counts Output File
        lastnt_outputFileName = os.path.join(outputDirectory, sampleName + "_" + sampleConc + "_LastntCounts.csv")
        if ReportNNNNNs == True:
            NNNNN_Tnc_outputFileName = os.path.join(outputDirectory, sampleName + "_" + sampleConc + "_NNNNN-Tnc.csv")
            NNNNN_FL_outputFileName = os.path.join(outputDirectory, sampleName + "_" + sampleConc + "_NNNNN-FL.csv")
            NNNNN_Combined_outputFileName = os.path.join(outputDirectory, sampleName + "_" + sampleConc + "_NNNNN-Combined.csv")
        
        # Determine behavior if output file already exists/open the outfile
        if os.path.isfile(lastnt_outputFileName):
            sys.stderr.write("Output file already exists; overwriting file!\n\n")
            Lastnt_fileOut = open(lastnt_outputFileName, "w")
            Lastnt_fileOut.write(",".join(str(val) for val in LastntHeaders) + '\n')
        else:
            Lastnt_fileOut = open(lastnt_outputFileName, "w")
            Lastnt_fileOut.write(",".join(str(val) for val in LastntHeaders) + '\n')
            
        
        if ReportNNNNNs == True:
            NNNNN_Tnc_fileOut = open(NNNNN_Tnc_outputFileName, "w")
            NNNNN_Tnc_fileOut.write(",".join(str(val) for val in NNNNN_Headers) + '\n')
            
            NNNNN_FL_fileOut = open(NNNNN_FL_outputFileName, "w")
            NNNNN_FL_fileOut.write(",".join(str(val) for val in NNNNN_Headers) + '\n')
            
            NNNNN_Combined_fileOut = open(NNNNN_Combined_outputFileName, "w")
            NNNNN_Combined_fileOut.write(",".join(str(val) for val in NNNNN_Headers) + '\n')
        
        
        # Setup main dictionary for the current file
        ReadsPFL_Dict[sampleConc] = {}
        
        # Add headers for the current file to header list
        headers.append(str(sampleConc) + '_Reads')
        headers.append(str(sampleConc) + '_PFL')
        
        # Adds molar conc. values for each file to conc_list
        conc_list.extend([0, molarConcValue])
            
        # Initialize some variables
        trunc_Dict = {}
        FL_Dict = {}
        early_Dict = {}
        header_counter = 0
        
        # Set file path
        variantPath = os.path.join(folder, file)
        
        # Read and analyze each specific file in the file list
        with open(variantPath, 'r') as f:
            # Generate events and quals for each SAM line
            for line in f:
                
                # Skip headers
                if header_counter == 0:
                    header_counter += 1
                    continue
                
                parsedLine = line.split('\t')
                #parsedLine = [ast.literal_eval(i) for i in parsedLine]            
                
                varName = ast.literal_eval(parsedLine[0])
                firstLast = ast.literal_eval(parsedLine[1])
                QNAME = parsedLine[2]
                if ReportNNNNNs == True:
                    NNNNN_val = parsedLine[3][:-1]
                    if 'N' in NNNNN_val:
                        NNNNN_val = 'XXXXX'
                
                lastntIndex = firstLast[1]-1 # making an index for use in the terminationDict (which is why 1 is being subtracted)
                
                sumReadsParsed += 1
                
                if len(varName) <= mutationsAllowed: # Only look at 0-[mutationsAllowed] mutations   
                    variantSet.add(varName)
                    
                    # Checking to see if the read spans the entire riboswitch 
                    # (otherwise it could get binned with the wrong variant)
                    sumReadsAnalyzed += 1
                    
                    
                    if firstLast[0] <= begin:
                        
                        if firstLast[1] >= trunc:
                            if varName[0][0] == 'WT':
                                zeroMutations += 1
                            elif len(varName) == 1:
                                oneMutation += 1
                            elif len(varName) == 2:
                                twoMutations += 1
                            else:
                                threeplusMutations += 1
                        # Keep track of the termination sites for looking at in a histogram
                        try:
                            LastntDict[varName][lastntIndex] += 1
                        except:
                            LastntDict[varName] = [0]*length
                            LastntDict[varName][lastntIndex] += 1
                            
                                
                        # Last nt is past the designated spot for FL reads            
                        if firstLast[1] >= FL:
                            try:
                                FL_Dict[varName] += 1
                            except:
                                FL_Dict[varName] = 1
                                trunc_Dict[varName] = 0
                                early_Dict[varName] = 0   
                            
                            if ReportNNNNNs == True:
                            # Keep track of NNNNN values for all 'full-length' RNAs
                                try:
                                    NNNNN_Combined_Dict[varName][NNNNN_indexDict[NNNNN_val]] += 1
                                    NNNNN_FL_Dict[varName][NNNNN_indexDict[NNNNN_val]] += 1
                                except:
                                    NNNNN_Combined_Dict[varName] = [0]*(Ns_counter+1)
                                    NNNNN_FL_Dict[varName] = [0]*(Ns_counter+1)
                                    NNNNN_Tnc_Dict[varName] = [0]*(Ns_counter+1)
                                    
                                    NNNNN_Combined_Dict[varName][NNNNN_indexDict[NNNNN_val]] += 1
                                    NNNNN_FL_Dict[varName][NNNNN_indexDict[NNNNN_val]] += 1
                                
                        # Last nt is between designated spots for Truncated and FL reads
                        elif firstLast[1] >= trunc:
                            try:
                                trunc_Dict[varName] += 1
                            except:
                                FL_Dict[varName] = 0
                                trunc_Dict[varName] = 1
                                early_Dict[varName] = 0
                            
                            # Keep track of NNNNN values for all 'truncated' RNAs
                            if ReportNNNNNs == True:
                                try:
                                    NNNNN_Combined_Dict[varName][NNNNN_indexDict[NNNNN_val]] += 1
                                    NNNNN_Tnc_Dict[varName][NNNNN_indexDict[NNNNN_val]] += 1
                                except:
                                    NNNNN_Combined_Dict[varName] = [0]*(Ns_counter+1)
                                    NNNNN_FL_Dict[varName] = [0]*(Ns_counter+1)
                                    NNNNN_Tnc_Dict[varName] = [0]*(Ns_counter+1)
                                    
                                    NNNNN_Combined_Dict[varName][NNNNN_indexDict[NNNNN_val]] += 1
                                    NNNNN_Tnc_Dict[varName][NNNNN_indexDict[NNNNN_val]] += 1
                                
                        # Last nucleotide comes before the termination site
                        # these reads are discarded because it is ambiguous if they should be classifed as FL or truncated
                        else:
                            sumEarlyTrunc += 1
                            try:
                                early_Dict[varName] += 1
                            except:
                                FL_Dict[varName] = 0
                                trunc_Dict[varName] = 0
                                early_Dict[varName] = 1
                    else:
                        sumBeganLate += 1
                        
                else:
                    if firstLast[0] <= begin and firstLast[1] >= trunc:
                        if len(varName) == 2:
                            twoMutations += 1
                        else:
                            threeplusMutations += 1
        
    
        sortedVariants = sorted(variantSet, key= lambda x:[(i[1],i[2]) for i in x if len(i) > 1])    
        sortedVariants = sorted(sortedVariants, key=len)    
    
        # Creates a dictionary where each key keeps track of the # reads and %FL for each variant
        for variant in sortedVariants:
            try:
                FL_reads = FL_Dict[variant]
            except:
                FL_reads = 0
        
            try:
                trunc_reads = trunc_Dict[variant]
            except:
                trunc_reads = 0
                
            total = FL_reads + trunc_reads
            if total != 0:
                PFL = 100*FL_reads/total
            else:
                PFL = 'NA'
            
            ReadsPFL_Dict[sampleConc][variant] = [total, PFL]
        
        
        mutTotal = 0
        percentZero = 0
        percentOne = 0
        percentTwo = 0
        percentThreePlus = 0
        mutTotal = zeroMutations + oneMutation + twoMutations + threeplusMutations
        if mutTotal != 0:
            percentZero = 100*zeroMutations/mutTotal
            percentOne = 100*oneMutation/mutTotal
            percentTwo = 100*twoMutations/mutTotal
            percentThreePlus = 100*threeplusMutations/mutTotal
        
        
        # Fix variant name and add it to the beginning of each line
        for variant in sortedVariants:  
            variantName = ''
            # Makes the variant names more readable
            for mutation in variant:
                for x in mutation:
                    if isinstance(x, int) == True: # Subtract the offset amount from the nt position
                        x -= offset
                        if x <= 0: # For DNA/RNA 0 isn't used as a positional value
                            x -= 1
                    variantName += str(x)
                if variant.index(mutation) != len(variant)-1: # Don't add ';' after last mutation
                    variantName += ';'
            line = [variantName]
            
            # Add last nt positions
            try:
                line.extend(LastntDict[variant])
            # When a given variant showed up elsewhere, but not for this particular concentration
            except:
                line.extend([0]*length)
            Lastnt_fileOut.write(",".join(str(val) for val in line) + '\n')
            
            if ReportNNNNNs == True:
                line_Ns_Tnc = [variantName]
                line_Ns_FL = [variantName]
                line_Ns_Combined = [variantName]
                
                # Add NNNNN values (Tnc)
                try:
                    line_Ns_Tnc.extend(NNNNN_Tnc_Dict[variant])
                # When a given variant showed up elsewhere, but not for this particular concentration
                except:
                    line_Ns_Tnc.extend([0]*(Ns_counter+1))
                NNNNN_Tnc_fileOut.write(",".join(str(val) for val in line_Ns_Tnc) + '\n')
                
                # Add NNNNN values (FL)
                try:
                    line_Ns_FL.extend(NNNNN_FL_Dict[variant])
                # When a given variant showed up elsewhere, but not for this particular concentration
                except:
                    line_Ns_FL.extend([0]*(Ns_counter+1))
                NNNNN_FL_fileOut.write(",".join(str(val) for val in line_Ns_FL) + '\n')
                
                # Add NNNNN values (Combined)
                try:
                    line_Ns_Combined.extend(NNNNN_Combined_Dict[variant])
                # When a given variant showed up elsewhere, but not for this particular concentration
                except:
                    line_Ns_Combined.extend([0]*(Ns_counter+1))
                NNNNN_Combined_fileOut.write(",".join(str(val) for val in line_Ns_Combined) + '\n')
            
        Lastnt_fileOut.close()
        if ReportNNNNNs == True:
            NNNNN_Tnc_fileOut.close()
            NNNNN_FL_fileOut.close()
            NNNNN_Combined_fileOut.close()
        
        
        print("Summary for file:", str(file))
        print("Reads parsed:", str(sumReadsParsed))
        print("Reads with an acceptable length:", mutTotal)
        print("\t0 Mutations:", str(zeroMutations) + ", ", str(round(percentZero,2)), "%")
        print("\t1 Mutation:", str(oneMutation) + ", ", str(round(percentOne,2)), "%")
        print("\t2 Mutations:", str(twoMutations) + ", ", str(round(percentTwo,2)), "%")
        print("\t3+ Mutations:", str(threeplusMutations) + ", ", str(round(percentThreePlus,2)), "%")
        print("Reads with fewer than %s mutation(s):" % (mutationsAllowed), sumReadsAnalyzed)
        print("\tReads discarded due to incorrect first nt position:", str(sumBeganLate))
        print("\tReads discarded due to early termination:", str(sumEarlyTrunc))
        print("Reads passing all filters:", (sumReadsAnalyzed-sumBeganLate-sumEarlyTrunc), "\n")        


    # Determine behavior if output file already exists/open the outfile
    if os.path.isfile(outputFileName):
        sys.stderr.write("Output file already exists; overwriting file!\n\n")
        fileOut = open(outputFileName, "w")
        fileOut.write(",".join(headers) + '\n')
        fileOut.write(",".join(str(val) for val in conc_list) + '\n')
    else:
        fileOut = open(outputFileName, "w")
        fileOut.write(",".join(headers) + '\n')
        fileOut.write(",".join(str(val) for val in conc_list) + '\n')
    
   
    for variant in sortedVariants:  
        variantName = ''
        # Makes the variant names more readable
        for mutation in variant:
            for x in mutation:
                if isinstance(x, int) == True: # Subtract the offset amount from the nt position
                    x -= offset
                    if x <=0: # For DNA/RNA 0 isn't used as a positional value
                        x -= 1
                variantName += str(x)
            if variant.index(mutation) != len(variant)-1: # Don't add ';' after last mutation
                variantName += ';'
        line = [variantName]
        
        # Add # of reads and percent Full Length for each concentration
        for conc in sampleConcList:
            try:
                line.extend(ReadsPFL_Dict[conc][variant])
            # When a given variant showed up elsewhere, but not for this particular concentration
            except:
                line.extend([0, 'NA'])
        fileOut.write(",".join(str(val) for val in line) + '\n') 
    fileOut.close()
    print("Complete! Output at %s and %s.\n" % (outputFileName, lastnt_outputFileName))
    
except:
    raise