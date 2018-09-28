
# Python script for creating heat maps based on the analysis performed using
# SMARTT-DAP (SMARTT - Data Analysis Pipeline)
#
# - Requires Python 3 (has been tested with version 3.5)
# - Public release 1.0
# - Copyright 2018 Chad Torgerson

###################################################################################
# GPL statement:
#
# This file is part of SMARTT - Heat Map Maker.
# 
# SMARTT - Heat Map Maker is free software: you can redistribute it  
# and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of 
# the License, or (at your option) any later version.
#
# SMARTT - Heat Map Maker is distributed in the hope that it will be 
# useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
###################################################################################

import re
import math
import operator
import functools
import os
import pandas as pd
import numpy as np
from bs4 import BeautifulSoup

try:   
    
    def rgb_to_hex(red, green, blue):
        if max(red, green, blue) > 255 or min(red, green, blue) < 0:
            raise ValueError('rgb value not in range(256)')
        else:
            return "#{0:02x}{1:02x}{2:02x}".format(red, green, blue)
    
    def apparEnergy(foldchange):
        ddG = R*T*math.log(foldchange)
        return ddG
    
    def determine_rgb_linear(rgb_color1, rgb_color2, Colorfract):
        interpolated_val = [0,0,0]    
        for i in range(3):
            interpolated_val[i] = (rgb_color1[i] * (1-Colorfract)) + (rgb_color2[i] * Colorfract)
        hex_color_code = rgb_to_hex(round(interpolated_val[0]), round(interpolated_val[1]), round(interpolated_val[2]))
        return hex_color_code
    
    # Determine arithmetic mean
    def mean(some_list):
        new_list = [0 if i < 0 else i for i in some_list]
        average = sum(new_list)/len(new_list)
        return average
    
    # Determine geometric mean
    def h_mean(some_list):
        new_list = [0 if i < 0 else i for i in some_list]
        recips_list = [1/x for x in new_list]
        harmonic_mean = len(recips_list)/sum(recips_list)
        return harmonic_mean
    
    # Determine geometric mean
    def g_mean(some_list):
        new_list = [0 if i < 0 else i for i in some_list]
        product_val = functools.reduce(operator.mul, new_list)
        geometric_mean = product_val ** (1/len(new_list))
        return geometric_mean
        
    def Linear_Change_Colorfract(SampleVal, slope, intercept):
        Colorfract = float((SampleVal-intercept)/slope)
        if Colorfract < 0:
            Colorfract = 0
        elif Colorfract > 1:
            Colorfract = 1
        return Colorfract     
    
    # Initialize variables
    R = 0.0019872036 # kcal K−1 mol−1     <-- determines dG units
    T = 310 # K    <-- temperature relevant for current experiment
    heat_range2 = None # Don't change this
    
    # Color spectrum options -- These values may be altered if desired
    rgb_black = (0, 0, 0)
    rgb_orange = (255, 102, 0)
    rgb_red = (255,26,26)
    rgb_white = (255,255,255)
    rgb_midgreen = (128, 204, 153)
    rgb_green = (0, 153, 51)
    rgb_purple = (224, 31, 144)
    rgb_yellow = (255, 204, 0)
    rgb_blue = (0, 51, 204)
    rgb_teal = (0, 255, 255)
    
    hex_black = '#000000'
    hex_teal = '#00FFFF'
    hex_blue = '#0066ff'
    hex_orange = '#ff6600'
    hex_red = '#FF1A1A'
    hex_grey = '#666666'
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~Variables to edit~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    # Indicating sample parameters
    # All options are currently case sensitive
    parameter_type = 'Ymin' # Options: 'K', 'Ymin', 'Ymax', 'Amp'
    value_type = 'mean' # Options: 'mean', 'best', 'worst', 'h_mean', 'g_mean' -- heatmap color will be based on the best, worst, or mean value of all mutants at each position
    heat_range = 'default' #Options: 'default', Some range of type: [5,55]    
    fits_generated_with = 'R'  # Options: 'PRSIM' or 'R'

    # File locations
    csvFolder = "demos"
    svgFolder = "demos"
    outFolder = "output"
    if fits_generated_with == 'PRISM':
        csvFile = os.path.join(csvFolder, "Cte_Fits_AllSingleMutations_PRISM.csv")
    elif fits_generated_with == 'R':
        csvFile = os.path.join(csvFolder, "Cte_Fits_AllSingleMutations_R.csv")
    svgFile = os.path.join(svgFolder, "Cte_annotated.svg")
    outFile = os.path.join(outFolder, "Cte_colored_" + parameter_type + '_' + value_type + ".svg")
    
    
    # Parameter identifiers (column/row titles) used in csv input file
    # ***NOTE: IF ONE OF THE PARAMETERS WASN'T USED WHEN GENERATING THE FITS, IT NEEDS 
    #          TO BE LABELLED AS 'None' BELOW***
    if fits_generated_with == 'PRISM':
        K = 'K'
        Ymin = 'Ymin'
        Ymax = 'Ymax'
        Amplitude = None
    elif fits_generated_with == 'R':
        K = 'K'
        Ymin = 'Ymin'
        Ymax = 'Ymax'
        Amplitude = 'Amp'
    
    # Setup Cutoffs for Min Amp for consideration in K heatplots
    minAmp_cutoff = 15 # Amp value needs to be greater than this for a K to be considered (~2x the Ymin of WT -- or as you judge appropriate)
    maxAmp_cutoff = 100 # Amp value needs to be less than this for a K to be considered

    # Setup heat map range for Colorfract function
    if heat_range == 'default':
        if parameter_type == 'Amp':    # WT = ~58
            heat_range = [10,58]    # Color 1 range
            heat_range2 = [58, 106] # Color 2 range
        elif parameter_type == 'K':    # WT = 1
            heat_range = [apparEnergy(25.7), apparEnergy(1)]    # Color 1 range
            heat_range2 = [apparEnergy(1), apparEnergy(1/25.7)] # Color 2 range
        elif parameter_type == 'Ymin':    # WT = ~7
            heat_range = [-43,7]    # Color 1 range
            heat_range2 = [7, 50]   # Color 2 range
        elif parameter_type == 'Ymax':    # WT = ~65
            heat_range = [15,65]    # Color 1 range
            heat_range2 = [65, 115] # Color 2 range
    
    if parameter_type == 'Amp' or parameter_type == 'Ymax' or parameter_type == 'K':
        rgb_low_color = rgb_red     # Color for values lower than WT
        rgb_mid_color = rgb_white   # Colow for WT
        rgb_high_color = rgb_blue   # Color for values higher than WT
    elif parameter_type == 'Ymin':
        rgb_low_color = rgb_red     # Color for values lower than WT
        rgb_mid_color = rgb_white   # Colow for WT
        rgb_high_color = rgb_blue   # Color for values higher than WT
        
        if heat_range2 == None:
            rgb_mid_color = rgb_white    # Colow for WT
            rgb_high_color = rgb_red   # Color for values higher than WT

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    
    # Setting up 
            
    # Decide maxORmin based on value_type designation
    if value_type == 'best':
        if parameter_type == 'Amp' or parameter_type == 'Ymax':
            maxORmin = 'max'
        elif parameter_type == 'K' or parameter_type == 'Ymin':
            maxORmin = 'min'
    elif value_type == 'worst':
        if parameter_type == 'Amp' or parameter_type == 'Ymax':
            maxORmin = 'min'
        elif parameter_type == 'K' or parameter_type == 'Ymin':
            maxORmin = 'max'
    elif value_type == 'mean' or value_type == 'h_mean' or value_type == 'g_mean':
        maxORmin = None
    else:
        print('value_type has unknown designation')

    # Determin slope and intercept for Colorfract function
    cf_slope = heat_range[1] - heat_range[0]
    cf_intercept = heat_range[0]
    
    if heat_range2 != None:
        cf_slope2 = heat_range2[1] - heat_range2[0]
        cf_intercept2 = heat_range2[0]
    
    # Initialize variables
    if fits_generated_with == 'PRISM':
        NUCS = ['A', 'U', 'C', 'G']
    elif fits_generated_with == 'R':
        NUCS = ['A', 'T', 'C', 'G']

    # Set regex patterns
    svgVariantsRegRE = re.compile(r"(^[A|U|T|C|G]{1})(-?[0-9]{1,3})(_dup[0-9]?)?")
    csvVaraintsRegRE = re.compile(r"(^[A|U|T|C|G]{1})(-?[0-9]{1,3})([A|U|T|C|G]{1})")
    legendRegRE = re.compile(r"(^color)([0-9]_)([0-9]{1,2}$)")
    legtextRegRE = re.compile(r"^legendT_([0-9]+)_([0-9]+)")
    legtitleRegRE = re.compile(r"^legend_title_([0-9]+)")
       
    # Set up style for zone outlines
    p_style = 'fill:#00ffff;fill-opacity:1;stroke:#000000;stroke-width:0.5;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;stroke-dasharray:none;stroke-dashoffset:0;display:inline'
    
    # Read csv file
    if fits_generated_with == 'PRISM':
        df=pd.read_csv(csvFile, sep=',', encoding = "ISO-8859-1", skipinitialspace=True, header=0, index_col = 0)
        df2=pd.DataFrame(df.values[2:5,:], columns=df.columns, index=df.index[2:5]).apply(pd.to_numeric)
    
    elif fits_generated_with == 'R':
        df=pd.read_csv(csvFile, sep=',', encoding = "ISO-8859-1", skipinitialspace=True, header=0, index_col = 0)
        df2=pd.DataFrame(df.values[:,0::2], columns=df.columns[0::2], index=df.index).apply(pd.to_numeric).T
    
    if Ymax == None:
        # Add row -- Maximum % FL (Ymax)
        Ymax = 'Ymax'
        Ymaxs = np.array(df2.loc[Ymin]) + np.array(df2.loc[Amplitude])
        Ymaxs = pd.DataFrame(Ymaxs).T
        Ymaxs.columns = df2.columns.values #
        df2.loc[Ymax,:] = Ymaxs.loc[0,:]
    elif Amplitude == None:
        # Add row -- Amplitude (Amp)
        Amplitude = 'Amp'
        Amps = np.array(df2.loc[Ymax]) - np.array(df2.loc[Ymin])
        Amps = pd.DataFrame(Amps).T
        Amps.columns = df2.columns.values #
        df2.loc[Amplitude,:] = Amps.loc[0,:]
    elif Ymin == None:
        # Add row -- Amplitude (Amp)
        Ymin = 'Ymin'
        Ymins = np.array(df2.loc[Ymax]) - np.array(df2.loc[Amplitude])
        Ymins = pd.DataFrame(Ymins).T
        Ymins.columns = df2.columns.values #
        df2.loc[Ymin,:] = Ymins.loc[0,:]

    # Add row -- Changes in apparent binding energy (ddG)
    FoldChanges = np.array(df2.loc[K])/np.array(df2['WT'][K])
    ddG_vfunction = np.vectorize(apparEnergy)
    ddGs = pd.DataFrame(ddG_vfunction(FoldChanges)).T
    ddGs.columns = df2.columns.values #
    df2.loc['ddG',:] = ddGs.loc[0,:]
        
    
    # Create Set of all Varaint Bases
    varSet = set()
    for variant in df2.columns.values:
        if variant == 'WT':
            continue
        else:
            try:
                csvVar = csvVaraintsRegRE.search(variant)
                csvVar2add = csvVar.group(1) + csvVar.group(2)
                varSet.update([csvVar2add])
            except:
                print('Create Dict & Set')
    
    
    # Parse xml with BeautifulSoup
    svg = open(svgFile, 'r').read()
    soup = BeautifulSoup(svg, features = 'xml')
    paths = soup.findAll('path')
    rects = soup.findAll('rect')
    texts = soup.findAll('text')
    
    for r in rects:
        if 'color' in r['id']:
            color_search = legendRegRE.search(r['id'])
            Colorfract = float(color_search.group(3))/10
            if color_search.group(2) == '1_':                
                fill_color = determine_rgb_linear(rgb_low_color, rgb_mid_color, Colorfract)
                r['style'] = re.sub(r'(fill:)#\w{6}', r'\1' + fill_color, p_style)
            elif color_search.group(2) == '2_':
                fill_color = determine_rgb_linear(rgb_mid_color, rgb_high_color, Colorfract)
                r['style'] = re.sub(r'(fill:)#\w{6}', r'\1' + fill_color, p_style)
        else:
            continue
    
    for t in texts:
        text_search = ""
        if 'text' in t['id']:
            continue
        if t['id'] == 'title':
            title_str = ' '.join([value_type.capitalize(), parameter_type, 'by Position'])
            t.contents[0].string = title_str
        elif 'legendT' in t['id']:
            text_search = legtextRegRE.search(t['id'])
            Colorfract = float(text_search.group(2))/10
            extreme1 = ''
            extreme2 = ''
            if text_search.group(1) == '1':
                if parameter_type == 'Amp' or parameter_type == 'Ymin' or parameter_type == 'Ymax':
                    units = '%'
                    if text_search.group(2) == '0':
                        extreme1 = 'â‰¤' # this weird symbol gets read as the 'less than or equal to' sign
                        extreme2 = 'â‰¥' # this weird symbol gets read as the 'greater than or equal to' sign
                else:
                    units = ''
                    if text_search.group(2) == '0':
                        extreme1 = 'â‰¥' # this weird symbol gets read as the 'greater than or equal to' sign
                        extreme2 = 'â‰¤' # this weird symbol gets read as the 'less than or equal to' sign
                t.contents[0].string = ' '.join([extreme1, str(round(cf_slope*Colorfract + cf_intercept,1)), units])
            
            elif text_search.group(1) == '2':
                if parameter_type == 'Amp' or parameter_type == 'Ymin' or parameter_type == 'Ymax':
                    units = '%'
                    if text_search.group(2) == '10':
                        extreme1 = 'â‰¤' # this weird symbol gets read as the 'less than or equal to' sign
                        extreme2 = 'â‰¥' # this weird symbol gets read as the 'greater than or equal to' sign
                else:
                    units = ''
                    if text_search.group(2) == '10':
                        extreme1 = 'â‰¥' # this weird symbol gets read as the 'greater than or equal to' sign
                        extreme2 = 'â‰¤' # this weird symbol gets read as the 'less than or equal to' sign
                t.contents[0].string = ' '.join([extreme2, str(round(cf_slope2*Colorfract + cf_intercept2,1)), units])


    # Get values for each variant
    for p in paths:
        SampleValList = []
        if ('path'in p['id']): # these are not named objects (ie. they are not associated w/ a nt position)
            continue
        elif('-' in p['id']): # also named objects (ie. they are not associated w/ a nt position)
            continue
        else:
            try:
                # Convert T to U in variant name (associated with object in svg)
                mutant = None
                mutant_Val = None
                
                mutSearch = svgVariantsRegRE.search(p['id'])
                
                if fits_generated_with == 'R':
                    if mutSearch.group(1) == 'U' or mutSearch.group(1) == 'u':
                        mutant_base = 'T' + mutSearch.group(2)
                    else:
                        mutant_base = mutSearch.group(1) + mutSearch.group(2)
                elif fits_generated_with == 'PRISM':
                        mutant_base = mutSearch.group(1) + mutSearch.group(2)
                
                # Create a list containing all 3 sample values
                for nt in NUCS:
                    mutant_Val = ''
                    if mutSearch.group(1) != nt:
                        mutant = mutant_base + nt
                        try:
                            if parameter_type == 'Amp':
                                mutant_Val = df2[mutant][Amplitude]
                            elif parameter_type == 'K':
                                if minAmp_cutoff < df2[mutant][Amplitude] < maxAmp_cutoff:
                                    mutant_Val = df2[mutant]['ddG']
                            elif parameter_type == 'Ymin':
                                mutant_Val = df2[mutant][Ymin]
                            elif parameter_type == 'Ymax':
                                mutant_Val = df2[mutant][Ymax]
                        except:
                            mutant_Val = ''
                        SampleValList.append(mutant_Val)
                
                # Remove all strings from SampleAmpList (such as "")
                SampleValList = [x for x in SampleValList if not isinstance(x, str)]
                if len(SampleValList) > 0:
                    if maxORmin == 'min':
                        SampleValue = min(SampleValList) # Determine largest Amp value at current nt position
                    elif maxORmin == 'max':
                        SampleValue = max(SampleValList)
                    elif value_type == 'mean':
                        SampleValue = mean(SampleValList)
                    elif value_type == 'h_mean':
                        SampleValue = h_mean(SampleValList)
                    elif value_type == 'g_mean':
                        SampleValue = g_mean(SampleValList)
                else:
                    SampleValue = ""
            except:
                print('whoops', mutant, mutant_Val, p['id'])
                continue
        
        # Set color for each position
        if SampleValue == '':           
            fill_color = hex_grey    
            p_style = re.sub(r'(fill:)#\w{6}', r'\1' + fill_color, p_style)
            p['style'] = p_style
            continue
        else: 
            # Rule for determining object color
            if parameter_type == 'Amp' or parameter_type == 'Ymax' or parameter_type == 'Ymin':
                if heat_range2 == None or SampleValue <= heat_range[1]:
                    Colorfract = Linear_Change_Colorfract(SampleValue, cf_slope, cf_intercept)
                    fill_color = determine_rgb_linear(rgb_low_color, rgb_mid_color, Colorfract)
                else:
                    Colorfract = Linear_Change_Colorfract(SampleValue, cf_slope2, cf_intercept2)
                    fill_color = determine_rgb_linear(rgb_mid_color, rgb_high_color, Colorfract)                   
            elif parameter_type == 'K':
                if heat_range2 == None or SampleValue <= heat_range[0]:
                    Colorfract = Linear_Change_Colorfract(SampleValue, cf_slope, cf_intercept)
                    fill_color = determine_rgb_linear(rgb_low_color, rgb_mid_color, Colorfract)
                else:
                    Colorfract = Linear_Change_Colorfract(SampleValue, cf_slope2, cf_intercept2)
                    fill_color = determine_rgb_linear(rgb_mid_color, rgb_high_color, Colorfract)            
            
            # Change fill colors
            p_style = re.sub(r'(fill:)#\w{6}', r'\1' + fill_color, p_style)
            p['style'] = p_style

        if mutant_base not in varSet:
            print('Warning:', p['id'], 'not found in csv file')

    # Output the edited SVG file
    f = open(outFile, "w")
    f.write(str(soup))
    f.close()
    
    print('Complete!')
    print('Parameter:\t', parameter_type)
    print('Value reported:\t', value_type.capitalize())
    
except:
    raise