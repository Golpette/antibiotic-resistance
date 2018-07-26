#################################################################################################
#
# Program to plot probability that a given drver mutation is present above some cutoff in each
# well. Idea is to be able to see the most probable mutational pathway (MPMP) and how this changes
# for different ciprofloxacin concentration profiles.
#
#################################################################################################


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



# IMPORTANT: CURRENTLY NEED TO SPECIFY THIS MANUALLY!! ######################
num_experiments = 200
#############################################################################



# Input directory containg the "well_x_...dat" files
path = join(".", str(sys.argv[1]))

# Define minimum abundance of PM to be experimentally detectable
CUTOFF = 0.1

# List of wells we will get stats for
#wells = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]
nwells = 24

# Define genomic position of driver PMs
# 0 = gyrA1,  1 = gyrA2,  2 = parC,  3 = marR,  4 = acrR
DRIVERS = [0,1,2,3,4]

# Store no. of occurrences when PM is present above cutoff in specified well
driv_0 = [0 for i in range(nwells)] # gyrA1
driv_1 = [0 for i in range(nwells)] # gyrA2
driv_2 = [0 for i in range(nwells)] # parC
driv_3 = [0 for i in range(nwells)] # marR
driv_4 = [0 for i in range(nwells)] # acrR




for w in range(nwells):

    for exp in range(0,num_experiments):


        # Root of files relating to current well
        search = "well_"+str(w)+"_(.*)_"+str(exp)+".dat"

        # Get file for given experiment and well 
        #required_files = [f for f in listdir('.') if( isfile(f) and re.search(search, f) ) ]
        required_files = [f for f in listdir(path) if( isfile(join(path,f)) and re.search(search, f) ) ] 

        if len(required_files) > 1: 
            print "ERROR: more than 1 file found for given experiment and well."
            print required_files, "###"

      
        for filename in required_files:
            
            f = open( join(path,filename), 'r' )
            for line in f:
                
                # Read in line as list of numbers (rank, PM_location, abundance) and store number
                # of experiments in which each driver has and abundance > CUTOFF
                columns = [ float(x) for x in line.split() ] 
                if ( columns[2] > CUTOFF and columns[1] in DRIVERS ):
                    
                    if( columns[1] == 0 ):
                        driv_0[w] += 1
                    if( columns[1] == 1 ):
                        driv_1[w] += 1
                    if( columns[1] == 2 ):
                        driv_2[w] += 1                         
                    if( columns[1] == 3 ):
                        driv_3[w] += 1
                    if( columns[1] == 4 ):
                        driv_4[w] += 1





######################  PLOTTING DATA FOR EACH DRIVER PM #######################################
# Plot probability each mutation is present in a given well, with each driver PM having a different shape/color

fig = plt.figure()
fig.set_size_inches(14, 8)
title = "Probability PM is present in each well (cutoff = "+str(CUTOFF)+")"
fig.suptitle(title, fontsize=32)
plt.xlabel('Well', fontsize=26);  
plt.ylabel('Probability', fontsize=22)
#plt.xlim((0,25));
#plt.ylim((0,1))
plt.tick_params(axis='both', which='major', labelsize=24)


# Make column scatter plot of "well vs #PMs" for all detectable PMs
w0=[];  p0=[];
w = 0
for n in driv_0:
  w0.append(w);   p0.append(float(n)/float(num_experiments));  w+=1;
lbl = "gyrA1"
plt.scatter(w0, p0, color='green', marker='s', s=500, label=lbl)

w1=[];  p1=[];
w = 0
for n in driv_1:
  w1.append(w);   p1.append(float(n)/float(num_experiments));  w+=1;
lbl = "gyrA2"
plt.scatter(w1, p1, color='red', marker='o', s=400, label=lbl)

w2=[];  p2=[];
w = 0
for n in driv_2:
  w2.append(w);   p2.append(float(n)/float(num_experiments));  w+=1;
lbl = "parC"
plt.scatter(w2, p2, color='black', marker='o', s=300, label=lbl)

w3=[];  p3=[];
w = 0
for n in driv_3:
  w3.append(w);   p3.append(float(n)/float(num_experiments));  w+=1;
lbl = "marR"
plt.scatter(w3, p3, color='blue', marker='o', s=200, label=lbl)

w4=[];  p4=[];
w = 0
for n in driv_4:
  w4.append(w);   p4.append(float(n)/float(num_experiments));  w+=1;
lbl = "acrR"
plt.scatter(w4, p4, color='orange', marker='o', s=100, label=lbl)



plt.legend(bbox_to_anchor=(1.03, 1), loc=2, borderaxespad=0.)

plt.show()
#figureName = "MPMP_cutoff_"+str(CUTOFF)+".png"
#plt.savefig( join(path,figureName) )
#plt.savefig( join(path,figureName), bbox_inches='tight')

################################################################################################







exit(0)

