# Program to:
# (1) make Fisher wave plots,
# (2) calculate speed of wavefront

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

# Plot Fisher wave for every "time_gap" number of data points
time_gap = 1



fig = plt.figure()
title = "Fisher waves"
fig.suptitle(title, fontsize=32)
plt.xlabel('Well', fontsize=26);  plt.ylabel('Well population', fontsize=26)
plt.tick_params(axis='both', which='major', labelsize=20)
plt.ylim((0,1.2))

#plt.yscale('log')


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


    # --- For plot of Fisher waves ---------------------------------
    line_num += 1
    if (line_num % time_gap == 0 ):

        time = columns[0]
        columns.pop(0)
        
        w = 0
        wells = []
        populations = []
        for i in columns:
            wells.append(w)
            w += 1
            populations.append( i )

        lbl = str(int(time))+" hours"
        plt.plot( wells, populations, linewidth='4', label=lbl )

    #---------------------------------------------------------------


plt.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)





######################  PLOTTING SPEED OF WAVEFRONT ############################################
fig2 = plt.figure()
title = "Position of wavefront vs time"
fig2.suptitle(title, fontsize=34)
plt.xlabel('t( N(well) > N* ) [hours]', fontsize=26);  plt.ylabel('Well', fontsize=26)
plt.tick_params(axis='both', which='major', labelsize=20)


# Make column scatter plot of "well vs #PMs" for all detectable PMs
w=[];  t=[];

abs_times.pop(0)  # Ignore first 2 data points or we see dependence on initial population size
abs_times.pop(0)

for (a,b) in abs_times:
    w.append(a);   t.append(b);

# Calculate slope of plot using numpy polyfit()
slope, intercept = np.polyfit( t, w, 1 )
speed = round_sig( slope, 3 )

text="Speed = "+str(speed)+" wells/h"
plt.plot(t, w, linewidth='3')
#plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
#plt.text(int(len(abs_times)*0.2), int(len(abs_times)*0.9), text, fontsize=26)
plt.text( (t[len(t)-1]*0.2), (w[len(w)-1]*0.9), text, fontsize=26)


plt.show()
################################################################################################








exit(0)

