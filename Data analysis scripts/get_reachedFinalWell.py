# Program to check how many runs reached final well by looking at "sequence files" which
# contain trajectory data for final well that had population > 0.1*K


# obtaining file lists
import os
from os import listdir
from os.path import isfile, join
# regexp; need search() method for substring in filename
import re
# plotting data
import numpy as np
import matplotlib.pyplot as plt
# arguments from command line
import sys


# Input directory containg the "well_x_...dat" files
path = join(".", str(sys.argv[1]))


# Define minimum abundance of PM to be experimentally detectable
CUTOFF = 0.1


# List of wells we will get stats for. #### MAKE SURE FINAL WELL IS FINAL ENTRY!! ####
#wells = [0, 3, 5, 7, 11, 13, 15, 19, 21, 23]
wells = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]

reached_final_well = 0
num_sims = 0


search2 = "sequences_\(well"
rq_files = [f for f in listdir(path) if( isfile(join(path,f)) and re.search(search2, f) ) ]
num_sims += len( rq_files )


final_well = wells[ len(wells)-1 ]
search = "sequences_\(well"+str(final_well)+"\)_"
# Make list of files relating to this well (ONLY from current directory '.')
required_files = [f for f in listdir(path) if( isfile(join(path,f)) and re.search(search, f) ) ]
reached_final_well = len(required_files)



print ""
print "Number runs = ", num_sims
print "Num that reached final well = ", reached_final_well



exit(0)

