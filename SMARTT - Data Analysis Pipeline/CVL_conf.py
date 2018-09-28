
# Config file for CreateVariantsList.py
# - Requires Python 3 (has been tested with version 3.5)
# - Public release 1.0
# - Copyright 2018 Chad Torgerson
#
# This file was adapted from the config file associated with the RTEventsCounter.py 
# script (DOI: 10.1021/acs.biochem.7b00323) by the Simon Lab at Yale. Thus, some
# of the options may be vestigial. 

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


# MAPQ threshold to include a SAM line
minMAPQ = 20

# Minimum quality score allowed for a given base. If any base has a score below
# this threshold, the line is skipped because the number/type of mutations is ambiguous.
qualityScoreThreshold = 20

# mutation Range -- for the C. tetani singlet: first mutation = 25; last = 198
mutRange = [25,198]


###################################################################################
###################################################################################
# Don't change any of the below options unless you know what you are doing.
# This code is vestigial from RTEventsCounter.py
###################################################################################
###################################################################################

# Overwrite or append when output file already exists
appendIfOutputFileExists = False

# Long deletions (>1nts)
# "dels"/"stops"/"none"
longDelsAs = "dels"

# Include insertions counts
countInserts = True

# Paired reads behavior
alignPaired = True
ignoreUnpaired = True

# Nucleotide numbering base in output (0-based or 1-based)
ntBase1 = True

# Output same-nucleotide toN as actual count or "NA"
sameNucAsNA = True
