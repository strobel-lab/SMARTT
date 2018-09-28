
# Python script for parsing through sequence alignments to create a list 
# of variants with their start and stop position for a library of mutant RNA 
# (Step 1 in the Pipeline)
#
# This code was adapted from the RTEventsCounter.py script 
# (DOI: 10.1021/acs.biochem.7b00323) by the Simon Lab at Yale. Thus, some of 
# the code is vestigial.
# 
# - Requires a config file (CVL_conf.py), a SAM file, and the reference genome 
#   used to create the SAM file
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


try:
    try:
        import CVL_conf
    except:
        sys.stderr.write("Error importing configuration file.\n")
        raise

    NUCS = ["A", "T", "G", "C"]
    COMP = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C"
    }

    def complement(seq):
        # Complement sequence
        return "".join([COMP.get(c,"N") for c in seq])

    def safeDivide(num, den):
        # Divide numbers as floats; -999 for zero division
        try:
            if num == "NA" or den == "NA":
                return "NA"
            else:
                return float(num)/float(den)
        except ZeroDivisionError:
            return -999.0
        except:
            sys.stderr.write("Unexpected error at division: " + str(sys.exc_info()[:2]) + "\n")
            return -999.0

    def checkSeq(line):
        # Check and fix FASTA sequence lines before parsing it
        uFound = False
        unknownFound = False
        for c in line:
            if c in ["a", "t", "g", "c"]:
                c = c.upper()
            elif c in ["U", "u"]:
                c = "T"
                uFound = True
            elif c not in NUCS:
                c = "N"
                unknownFound = True
        if uFound:
            sys.stderr.write("FASTA notice: U found and substituted as T\n")
        if unknownFound:
            sys.stderr.write("FASTA warning: Non-nucleotide character found!\n")
        return line

    def parseFasta(fastaFile):
        # Takes a FASTA file and returns a dict of sequences associated with their names
        seqs = {}
        seqName = ""
        with open(fastaFile, 'r') as f_fasta:
            for line in f_fasta:
                if line[0] == ">":
                    # Sequence name
                    seqName = line.strip()[1:]
                    if seqName in seqs:
                        sys.stderr.write("FASTA error: Duplicate seq name ("+seqName+") found in FASTA!\n")
                    else:
                        seqs[seqName] = ""
                elif seqName != "":
                    # Sequence, strip whitespace
                    seqs[seqName] += checkSeq("".join(line.split()))
                    seqName = ""
                else:
                    sys.stderr.write("FASTA error: Skipped FASTA line, no seq name associated.\n")
        return seqs

    def parseCigarString(cigarString, read, quals, refSeq, startIndex):
        # Parse CIGAR string and generate raw events and quals
        # Resolve CIGAR strings into lists of numbers and event types
        cigarRegRE = re.compile(r"[0-9]+[M|I|D|S]")  # Note: Only CIGAR strings with M, I, D, S are allowed
        numberRE = re.compile(r"(^[0-9]+)")
        splitCigarString = cigarRegRE.findall(cigarString.strip())
        cigarList = [numberRE.split(reg)[1:] for reg in splitCigarString]

        # Prepare output dicts
        events = {}
        alignedQuals = {}

        # Prepare variables
        lastQual = "#"  # For assigning qual to deleted region
        refIndex = startIndex  # Begin at leftmost matching position; 0-based

        # Go through every CIGAR region
        try:
            for i in range(len(cigarList)):
                region = cigarList[i]
                regionLength = int(region[0])
                regionType = region[1]

                # Match/mismatch
                if regionType == "M":
                    matchRegion = read[:regionLength]
                    qualsRegion = quals[:regionLength]
                    read = read[regionLength:]
                    quals = quals[regionLength:]
                    for regionIndex in range(regionLength):
                        nuc = matchRegion[regionIndex]
                        qual = qualsRegion[regionIndex]
                        lastQual = qual
                        alignedQuals[refIndex] = qual
                        if nuc != refSeq[refIndex]:
                            events[refIndex] = nuc
                        else:
                            events[refIndex] = "|"
                        refIndex += 1
                # Insertion
                elif regionType == "I":
                    if (CVL_conf.countInserts and
                            events[refIndex-1] == "|" and refSeq[refIndex] == read[regionLength] and
                            read[0] != refSeq[refIndex] and read[regionLength-1] != refSeq[refIndex-1]):
                        # Record on 5 prime neighboring nuc if unambiguous, and not flanked by mismatches
                        events[refIndex-1] = "^"
                    read = read[regionLength:]
                    quals = quals[regionLength:]
                # Deletion
                elif regionType == "D":
                    for delIndex in range(refIndex, refIndex+regionLength):
                        events[delIndex] = "-"
                        alignedQuals[delIndex] = lastQual
                    refIndex += regionLength
                # Soft clipping
                elif regionType == "S":
                    try:
                        read = read[regionLength:]
                        quals = quals[regionLength:]
                    except IndexError:
                        pass
                    clipRange = []
                    if i == len(cigarList)-1:  # Rightmost end of read/CIGAR
                        clipRange = range(refIndex, refIndex+regionLength)
                    elif i == 0:  # Leftmost end of read/CIGAR
                        clipRange = range(refIndex, refIndex-regionLength-1, -1)
                    if clipRange != []:
                        for clipIndex in clipRange:
                            if clipIndex >= 0 and clipIndex < len(refSeq):
                                events[clipIndex] = "s"
                                alignedQuals[clipIndex] = "#"
        except (IndexError, KeyError):
            sys.stderr.write("CIGAR parsing index/key error!\n")
            return None, None

        return events, alignedQuals

    def cullEvents(events, refSeq):
        # cull Events
        eventChars = ["s", "-", "A", "T", "G", "C", "N", "^"]
        culledEvents = dict(events)        
        
        # Scan left to right
        for i in range(min(events.keys()), max(events.keys())+1):
            if i not in events.keys():
                # Absent gaps in events
                culledEvents[i] = "~"
            
            # This section used to prevent consecutive mutations, but they should be allowed here
            # Not sure why I can't just remove the section completely, but if I do then no mutations are written
            # To the VariantList file
            elif events[i] in NUCS:
                # Mismatch
                try:
                    if events[i+1] in eventChars:
                        pass
#                        # With adjacent downstream mutation events
#                        culledEvents[i] = "|"
                except (KeyError, IndexError):
                    pass
            elif events[i] == "-":
                # Deletion
                if events[i+1] not in eventChars:
                    # With no adjacent downstream mutation events (|/~)
                    # Handle ambiguous deletions
                    k = 1
                    while events[i-k] == "-":
                        k += 1
                    fivePrimeStart = i-k+1
                    threePrimeEnd = i
                    ambiguous = False
                    # Slide deletion upstream and downstream
                    notDelSeq = refSeq[:fivePrimeStart] + refSeq[threePrimeEnd+1:]
                    maxOffset = 1 + threePrimeEnd-fivePrimeStart
                    for offset in range(-maxOffset, maxOffset+1):
                        notSubSeq = ""
                        if offset != 0:
                            try:
                                notSubSeq = refSeq[:fivePrimeStart + offset] + refSeq[threePrimeEnd + offset+1:]
                            except (KeyError, IndexError):
                                pass
                        if notSubSeq == notDelSeq:
                            ambiguous = True
                    if not ambiguous:
                        culledEvents[i] = "-"

                elif events[i+1] != "-":
                    # At the edge with soft clip (s), or with downstream mutation (A/T/G/C/N), or insert (^)
                    culledEvents[i] = "~"
                else:
                    # Long deletion (-)
                    culledEvents[i] = "~"
                    if events[i+2] != "-":
                        # start of long deletion at i+1; record as stop if opted
                        if CVL_conf.longDelsAs == "stops":
                            culledEvents[i+1] = "F"
                        elif CVL_conf.longDelsAs == "none":
                            culledEvents[i+1] = "~"
            #elif events[i] == "^" or "~" or "s" #insertion or __ or softclip
            elif events[i] not in ["~", "|", "s", "^"]:
                # Unwanted events; replace with match character
                culledEvents[i] = "|"   
                        
        # Record start/stop events
        firstNt = min(culledEvents.keys())
        # c= start site
        lastNt = max(culledEvents.keys())
        # c= end site
        
        # Skip over soft clips
        if culledEvents[firstNt] == "s":
            while culledEvents[firstNt] == "s":
                firstNt += 1        
        
        # Skip over soft clips
        if culledEvents[lastNt] == "s":
            while culledEvents[lastNt] == "s":
                del culledEvents[lastNt]
                lastNt -= 1
                if len(culledEvents) == 0:
                    break
        if lastNt < len(refSeq)-1:
            culledEvents[lastNt+1] = "!"
               
        firstLastNt = (firstNt+1,lastNt+1)        
        return culledEvents, firstLastNt

    def nameVariant(events, refSeq):
        # nameVariant
        numMuts = 0
        mutName = (('WT',),)
        mutsOnly = True     
        mut_lower = CVL_conf.mutRange[0]
        mut_upper = CVL_conf.mutRange[1]
        
        # Scan left to right
        for i in range(min(events.keys()), max(events.keys())+1):
            if i not in events.keys():
                # Absent gaps in events
                #print(i)                
                mutsOnly = False
            elif events[i] in NUCS:
                # Mismatch
                tup = (refSeq[i], int(i+1), events[i])
                
                if i >= (mut_lower-1) and i < mut_upper:                
                    try:
                        mutName.append(tup)
                    except:
                        mutName = [tup]
                    numMuts += 1
            elif events[i] not in ["|", "!"]:
                # Deletion, long deletion, soft clip, insertion
                mutsOnly = False

        mutName = tuple(mutName)
        return mutName, numMuts, mutsOnly

    def checkQualityScores(events, qualities, threshold):
        # check to see if all quality scores are above threshold
        #count = 0
        AboveThreshold = True
        mut_lower = CVL_conf.mutRange[0]
        mut_upper = CVL_conf.mutRange[1]
        
        # Scan left to right
        for i in range(min(events.keys()), max(events.keys())):
            if i >= (mut_lower-1) and i < mut_upper:
                try:
                    score = ord(qualities[i])-33
                    if score < threshold:
                        AboveThreshold = False
                except:
                    AboveThreshold = False

        return AboveThreshold


    def parseSamLine(line):
        # Takes a SAM file line of alignment and returns a list of the components
        splitLine = line.split()
        # Get relevant fields
        try:
            refName = splitLine[2]
            if refName not in refSeqs.keys():
                return ["skip", None, None, None, None, None]
            clusterName = splitLine[0]
            startIndex = int(splitLine[3])-1  # 0-based
            mapQual = int(splitLine[4])
            cigarString = splitLine[5]
            rawRead = splitLine[9]
            rawQual = splitLine[10]
        except IndexError:
            sys.stderr.write("SAM line parsing index error.\n")
            return ["skip", None, None, None, None, None]

        if cigarString != "*":
            events, qualities = parseCigarString(cigarString, rawRead, rawQual, refSeqs[refName], startIndex)
            if events is not None:
                return [clusterName, refName, events, qualities, mapQual, rawRead]
            else:
                return ["skip", None, None, None, None, None]
        else:
            return ["skip", None, None, None, None, None]

    def combineEvents(R1, R2, Q1, Q2):
        # Combine two events/quals dicts for F/R seqs
        combinedEvents = dict(R1)
        combinedQuals = dict(Q1)
        # overlappingNucs = []
        for i in R2.keys():
            combinedEvents[i] = R2[i]
            combinedQuals[i] = Q2[i]
            if i in R1.keys():
                # overlappingNucs.append(i)
                combinedQuals[i] = chr(int((ord(Q1[i])+ord(Q2[i]))/2.0))
                # Handle positions that do not agree between mate pairs
                if R1[i] != R2[i]:
                    try:
                        # Disregard soft clips
                        if R1[i] == "s":
                            combinedEvents[i] = R2[i]
                        elif R2[i] == "s":
                            combinedEvents[i] = R1[i]
                        # Take higher qual events
                        elif ord(Q2[i])-33 > ord(Q1[i])-33:
                            combinedEvents[i] = R2[i]
                        else:
                            combinedEvents[i] = R1[i]
                    except IndexError:
                        print("len(R1) = %i, len(Q1) = %i, len(R2) = %i, len(Q2) = %i" % (len(R1), len(Q1), len(R2), len(Q2)))
                        raise
        return combinedEvents, combinedQuals
    
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~Variables that can be edited~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # If debugging, set to True and file inputs will be ignored
    IgnoreFileInputs = True
    
    # These must be set if IgnoreFileInputs is set to True
    # Files and input data
    if IgnoreFileInputs == False:
        demoDirectory = "demos"
        samFilePath = os.path.join(demoDirectory, "Cte_100mM_Gly-demo.sam")
        refFilePath = os.path.join(demoDirectory, "Cte_ref.fa")
        outDirectory = "output"
    else:
        outDirectory = "output"
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
    
    
    # Check file paths inputs
    if len(sys.argv) != 3:
        InputLength = False
    else:
        InputLength = True
    
    if InputLength == False and IgnoreFileInputs == False:
        sys.exit("Usage: python MakeVariantList.py <sam file> <reference fasta>")
    else:          
        # Files and input data
        if IgnoreFileInputs == False:
            samFilePath = sys.argv[1]
            refFilePath = sys.argv[2]
        
        # Save sample name; remove path and filename extension
        sampleName = os.path.splitext(os.path.split(samFilePath)[1])[0]

        # Save output file names
        outputFileName_variants = os.path.join(outDirectory, sampleName + "_variants.txt")

        # Get reference seqs
        refSeqs = parseFasta(refFilePath)

        # Prepare output                
        header_variants = "\t".join([
            "variant",
            "(start, stop)",
            "QNAME"+ "\n"
            ])
 
        # Determine behavior if output file2 already exists
        if os.path.isfile(outputFileName_variants):
            if CVL_conf.appendIfOutputFileExists:
                sys.stderr.write("Output file already exists; appending current output to file!\n\n")
                fileOut_variants = open(outputFileName_variants, "a")
            else:
                sys.stderr.write("Output file already exists; overwriting file!\n\n")
                fileOut_variants = open(outputFileName_variants, "w")
                fileOut_variants.write(header_variants)
        else:
            fileOut_variants = open(outputFileName_variants, "w")
            fileOut_variants.write(header_variants)
           

        # Prepare data dicts
        mismatch = {}  # Mismatch event: mismatch[seqName][fromNuc][toNuc][i]
        mutationCount = {}  # Total mut count across mismatch types: mutationCount[seqName][i]
        deletion = {}  # Deletion events: deletion[seqName][fromNuc][i]
        stops = {}  # Stop events: stops[seqName][i]
        insertion = {}  # Insertion events: insertion[seqName][i]
        depth = {}  # Coverage depth: depth[seqName][i]
        RT = {}  # RT read through: RT[seqName][i]
        qualsList = {}  # Average alignment quals per nuc: qualsList[seqName][i][...]


        # For summary at the end
        sumSamLines = 0
        sumSamLowQual = 0
        sumSamDiscard = 0
        sumUnpairedReads = 0
        sumAboveQualThreshold = {}
        sumBelowQualThreshold = {}
        sumIndelsObserved = {}
        sumCulledEventsPerGroup = {}
        sumPassedFilters = {}
        for seqName in refSeqs:
            sumCulledEventsPerGroup[seqName] = 0
            sumAboveQualThreshold[seqName] = 0
            sumBelowQualThreshold[seqName] = 0
            sumIndelsObserved[seqName] = 0
            sumPassedFilters[seqName] = 0

        # Build data structure
        for seqName in refSeqs:
            mismatch[seqName] = {}
            for fromNuc in NUCS:
                mismatch[seqName][fromNuc] = {}
                for toNuc in [mutNuc for mutNuc in NUCS if mutNuc != fromNuc]:  # Only dict for mismatches
                    mismatch[seqName][fromNuc][toNuc] = [0]*len(refSeqs[seqName])
            mutationCount[seqName] = [0]*len(refSeqs[seqName])
            deletion[seqName] = {}
            for fromNuc in NUCS:
                deletion[seqName][fromNuc] = [0]*len(refSeqs[seqName])
            stops[seqName] = [0]*len(refSeqs[seqName])
            insertion[seqName] = [0]*len(refSeqs[seqName])
            depth[seqName] = [0]*len(refSeqs[seqName])
            RT[seqName] = [0]*len(refSeqs[seqName])
            qualsList[seqName] =  [[] for n in range(len(refSeqs[seqName]))]
            
        # Go through sequences and produce counts
        covered = ["|", "-", "A", "T", "G", "C", "^"]  # for depth
        readThrough = ["|", "A", "T", "G", "C", "^", "~"]  # for RT
        # For recording previous data for pair combination
        prevClusterName = ""
        prevTargetName = ""
        prevLine = ""
        prevEvents = []
        prevQualities = []
        prevRawRead = ""
        # Control booleans
        pairFound = False
        updateCounts = False        

        
        # Print sample name for summary
        print("Sample:", sampleName)   
       
        # Declare Refseq names
        print("Reference sequence names:")
        for seqName in refSeqs:
            print("\t" + str(seqName))
        
        print("\nParsing and analyzing events...")
        
        # Read sam files
        with open(samFilePath, 'r') as f:
            # Generate events and quals for each SAM line
            for line in f:
                if line[0] == "@":
                    # Skip headers
                    continue
                parsedLine = parseSamLine(line)
                sumSamLines += 1
                clusterName = parsedLine[0]
                targetName = parsedLine[1]
                events = parsedLine[2]
                qualities = parsedLine[3]
                mappingQual = parsedLine[4]
                rawRead = parsedLine[5]
    
                isCombinedRead = False
                shouldCount = False
                
                
    
                if clusterName == "skip":
                    # Loop through each line of data
                    sumSamDiscard += 1
                elif mappingQual < CVL_conf.minMAPQ:
                    # Exclude reads whose mapping quals are too low
                    sumSamLowQual += 1
                else:
                    # First, handle events and combine pair reads
                    if CVL_conf.alignPaired:
                        if not pairFound:
                            if clusterName == prevClusterName and targetName == prevTargetName:
                                # Found matching pair; combine reads
                                pairFound = True
                                updateCounts = True
                                isCombinedRead = True
                                
                                eventsToWrite, qualsToWrite = combineEvents(prevEvents, events, prevQualities, qualities)
                                # sortedKeys = sorted(eventsToWrite.keys())
                                
                            elif prevClusterName != "":
                                # Non-matching read
                                # Update with prev info if counting unpaired, but not in the 1st line case (no 0th line)
                                sumUnpairedReads += 1
                                if not CVL_conf.ignoreUnpaired:
                                    # overlappingNucs = []
                                    eventsToWrite = dict(prevEvents)
                                    qualsToWrite = dict(prevQualities)
                                    pairFound = False
                                    updateCounts = True
                        else:
                            # Pair found
                            pairFound = False
                            updateCounts = False
                    elif prevClusterName != "" and not CVL_conf.ignoreUnpaired:
                        updateCounts = True
                        eventsToWrite = dict(prevEvents)
                        qualsToWrite = dict(prevQualities)
                        # overlappingNucs = []
    
                    # Second, cull events dicts and update counts
                    if updateCounts:
                        updateCounts = False
    
                        # Cull events
                        refSeq = refSeqs[prevTargetName]
                        culledEventsToWrite, firstLast = cullEvents(eventsToWrite, refSeq)
                        
                        # Check to see if quality score for all positions is above the threshold (Q score of 20 = 99% accuracy)
                        goodQuality = False
                        goodQuality = checkQualityScores(culledEventsToWrite, qualsToWrite, CVL_conf.qualityScoreThreshold)
                        
                        if goodQuality == True:
                            sumAboveQualThreshold[prevTargetName] += 1
                            seqNameToCount, numMutations, shouldCount = nameVariant(culledEventsToWrite, refSeq)
                            
                            # Check to see if only mutations (not insertions or deletions) are present
                            if shouldCount == True: 
                                if numMutations >= 0:
                                    sumPassedFilters[prevTargetName] += 1
                                    
                                    fileOut_variants.write("\t".join([
                                        str(seqNameToCount),
                                        str(firstLast),
                                        str(clusterName)+ "\n"
                                        ]))
        
                                else:
                                    print('Warning: Negative number of mutations observed')
                            else:
                                sumIndelsObserved[prevTargetName] += 1
                        else:
                            sumBelowQualThreshold[prevTargetName] += 1
                        
                        # Update number of reads per group
                        sumCulledEventsPerGroup[prevTargetName] += 1
                        
                        
    
    
                    # Record data for pairing
                    prevClusterName = clusterName
                    prevTargetName = targetName
                    prevLine = line
                    prevEvents = events
                    prevRawRead = rawRead
                    prevQualities = qualities



        print("Complete! Output at %s.\n" % outputFileName_variants)

        # Print summary information
        print("Summary of run:")
        print("Each half of a paired-end reads counted separately:")
        print("\t" + "Reads parsed:", sumSamLines)
        print("\t" + "MAPQ score threshold:", CVL_conf.minMAPQ)
        print("\t" + "Reads skipped due to low MAPQ score:", sumSamLowQual)
        print("\t" + "Reads discarded (no ref. sequence, missing SAM fields, etc):", sumSamDiscard)
        print("\t" + "Unpaired reads:", sumUnpairedReads)
        
        print("\n" + "Paired-end reads combined before counting:")
        for seqName in refSeqs:
            print("\t" + "Reference Sequence:", seqName)
            print("\t" + "Reads passing first round of filters:", sumCulledEventsPerGroup[seqName])
            print("\t" + "Quality score threshold (QST) for individual bases:", CVL_conf.qualityScoreThreshold)
            print("\t" + "Quality reads (scores of all bases in mut. range > QST):", sumAboveQualThreshold[seqName])
            print("\t" + "Quality reads discarded for containing indels:", sumIndelsObserved[seqName])
            print("\t" + "Reads passing all filters:", sumPassedFilters[seqName], "\n")
        
        fileOut_variants.close()
except:
    raise
