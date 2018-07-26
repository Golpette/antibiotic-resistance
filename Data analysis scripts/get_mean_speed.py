# Program to calculate average speed of wavefronts
# across multiple experiments.
#
# NOTE: all "runs", in nwell.dat, must start at time=0
#


import numpy as np
# arguments from command line
import sys
# for round_sig() method
from math import log10, floor


def round_sig(x, sig):
    ''' Method to round to specified sig figs '''
    return round(x, sig-int(floor(log10(x)))-1)



# Arbitrary well occupancy for mid-point of wave 
n_star = 50000




filename = str(sys.argv[1])
f = open( filename, 'r' )


# Store absolute time intervals between colonisation of wells
absolute_ints = []
# and the absolute times
abs_times = []

well = 1
wave_front = 1
line_number = 0

speeds = []

for line in f:

    line_number+=1
    # Read in line as list of floats
    columns = [ float(x) for x in line.split()]
    # Number of wells (first column=time, all rest are wells)
    nwells = len(columns)-1


    # If time=0, new experiment, store speed and reset ---------
    if ( columns[0]==0 and line_number>1 ):

        #  Avg speed from fitting line to data ###################################
        w=[];  t=[];
        for (a,b) in abs_times:
            w.append(a);   t.append(b);

        # Calculate slope of plot using numpy polyfit()
        slope, intercept = np.polyfit( t, w, 1 )
        speed = round_sig( slope, 3 )

        speeds.append( speed )

        # Reset counters
        absolute_ints[:] = [] 
        abs_times[:] = []     
        wave_front = 1  
    # ----------------------------------------------------------


  
    # For generating gradient plot of speed --------------------
    if (wave_front < nwells and columns[ wave_front ] > n_star ):
         abs_times.append(   (wave_front, columns[0])   )
         wave_front += 1
    # ----------------------------------------------------------



# Rememebr to add speed from final data set
w=[];  t=[];

# Ignore first 2 points 
abs_times.pop(0)
abs_times.pop(0)

for (a,b) in abs_times:
    w.append(a);   t.append(b);
# Calculate slope of plot using numpy polyfit()
slope, intercept = np.polyfit( t, w, 1 )
speed = round_sig( slope, 3 )
speeds.append( speed )





print speeds
# Average time between wells
avg_speed = round_sig( np.mean( speeds ), 3 )   
std_speed = round_sig( np.std( speeds ), 3 )
print "Mean speed = ", avg_speed, "(+/-", std_speed,") wells/hour"
print ""





exit(0)

