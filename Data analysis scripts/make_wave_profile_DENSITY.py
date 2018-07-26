# Program to make profile plots of OD vs z=x-vt like Bartek does
#
# Note that it does this in real space (x in mm, not wells) and
# thus calculated speed in wells/h is transformed to mm/hour too
#

import numpy as np
import matplotlib.pyplot as plt
# arguments from command line
import sys
# for round_sig() method
from math import log10, floor



def round_sig(x, sig):
    ''' Method to round to specified sig figs '''
    return round(x, sig-int(floor(log10(x)))-1)



# Arbitrary well occupancy for mid-point of wave 
n_star = 0.5

# Get wave profile for every "time_gap" number of data points
time_gap = 50

# Need these to convert to real-space coords
#well_width = 3.25
#well_sep = 1.0






# FIRST: GET WAVESPEED FROM nwell.dat FILE ============================================
filename = str(sys.argv[1])
f = open( filename, 'r' )

# Store absolute time intervals between colonisation of wells
absolute_ints = []
# and the absolute times
abs_times = []

well = 1
wave_front = 1
line_num = 0

for line in f:

    # Read in line as list of floats
    columns = [ float(x) for x in line.split()]
    # Number of wells (first column=time, all rest are wells)
    nwells = len(columns)-1
    
    # For generating gradient plot of speed --------------------
    if (wave_front < nwells and columns[ wave_front ] > n_star ):
         abs_times.append(   (wave_front, columns[0])   )
         wave_front += 1
    # ----------------------------------------------------------

# --- Get wavespeed as gradient of wavefront vs time
# Make column scatter plot of "well vs #PMs" for all detectable PMs
w=[];  t=[];
for (a,b) in abs_times:
    w.append(a);   t.append(b);

# Calculate slope of plot using numpy polyfit()
slope, intercept = np.polyfit( t, w, 1)
speed = round_sig( slope, 3 )

print "speed = ", speed, "wells/hour"

# ======================================================================================






# SECOND: Use speed to make population vs z=x-vt plot =================================

filename = str(sys.argv[1])
f = open( filename, 'r' )

# Store absolute time intervals between colonisation of wells
absolute_ints = []
# and the absolute times
abs_times = []

line_num = 0

for line in f:

    # Read in line as list of floats
    columns = [ float(x) for x in line.split()]
    # Number of wells (first column=time, all rest are wells)
    nwells = len(columns)-1


    # -Collapse Fisher waves (ony if final well empty)----------------
    line_num += 1
    if (line_num % time_gap == 0   and  columns[nwells]<0.3):

        time = columns[0]
        columns.pop(0)
        
        w = 0
        wells = []
        populations = []
        for i in columns:
            wells.append(w)
            w += 1
            populations.append( i )

        # Print profile as population vs z=x-vt in REAL SPACE.
        for w in wells:
            
            x_coord = w
            #x_coord = (well_width/2.0) + ( w * (well_width+well_sep) )  ### mid-point + well_num*(width+chanel_length)
            real_speed = speed
            #real_speed = speed * (well_width+well_sep)

            z_coord = x_coord - (real_speed * time) 
            print z_coord, " ", populations[w]
    #---------------------------------------------------------------






exit(0)

