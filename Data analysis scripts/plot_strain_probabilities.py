#################################################################################################
#
# Program to plot probability that a given MUTANT STRAIN is present above some cutoff in each
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
# for colours
import matplotlib.cm as cm



# IMPORTANT: CURRENTLY NEED TO SPECIFY THIS MANUALLY!! ######################
num_experiments = 100

num_strains = 32
#############################################################################



## QUICK FIX: I was naming strains with integers before. Stupid. Want to print actual genotype as label
strainMap = {0:'00000', 1:'10000',2:'01000',3:'00100',4:'00010',5:'00001',
             6:'11000',7:'10100',8:'10010',9:'10001',10:'01100',11:'01010',12:'01001',13:'00011',
             14:'11100',15:'11010',16:'11001',17:'10110',18:'10101',19:'01110',20:'01101',21:'10011',22:'01011',
             23:'11011',24:'11110',25:'11101',26:'10111',27:'01111',28:'11111',
             29:'00110',30:'00101',31:'00111'}
#### This is hardcoded in the main program and is the order they appear in Marcusson2009 (+3 extras at end)





# Input directory containing the "strain_wx___...dat" files
path = join(".", str(sys.argv[1]))

# Define minimum abundance of strain to be experimentally detectable  SHOULD THIS BE n/K or n/N?? i.e. later wells might have very small population, <10% of K 
CUTOFF = 0.1
#CUTOFF=500

nwells = 24




# Store no. of occurrences of each strain in each well
strain_abunds = []
for n in range(num_strains):
    strain_abunds.append( [0 for i in range(nwells)]    )



for w in range(nwells):

    for exp in range(0,num_experiments):


        # Root of files relating to current well
        search = "strains_w"+str(w)+"___(.*)_"+str(exp)+".dat"

        # Get file for given experiment and well 
        required_files = [f for f in listdir(path) if( isfile(join(path,f)) and re.search(search, f) ) ] 

        if len(required_files) > 1: 
            print "ERROR: more than 1 file found for given experiment and well."
            print required_files, "###"

      
        for filename in required_files:
            
            f = open( join(path,filename), 'r' )
            for line in f:
                
                # Read in line as list of numbers (strain, abundance) and store number
                # of experiments in which each strain has an abundance > CUTOFF
                columns = [ float(x) for x in line.split() ] 

                strain = int(columns[0])

                if ( columns[1] > CUTOFF ):
                    # Increment counter
                    strain_abunds[strain][w] += 1
                    





#######################################################################################
# Plot probability each STRAIN is present in a given well, above cutoff abundance

fig = plt.figure()
fig.set_size_inches(14, 8)
title = "Probability strain is present in well (cutoff = "+str(CUTOFF)+")"
fig.suptitle(title, fontsize=32)
plt.xlabel('Well', fontsize=26);  
plt.ylabel('Probability', fontsize=26)
plt.xlim((0,25));
plt.ylim((0,1.1))
plt.tick_params(axis='both', which='major', labelsize=24)

# Iterate through color space
colors = iter(   cm.rainbow( np.linspace(0, 1, num_strains) )   )
#colors = iter(   cm.rainbow( np.linspace(0, 1, 15) )   )

for s in range(num_strains):
    strain_list = strain_abunds[s]
    w=[];  p=[];  
    #lbl="Strain "+str(s);
    lbl=strainMap[s]
    for wn in range(nwells):
        w.append( wn )
        p.append( float(strain_list[wn]) / float(num_experiments)   )

    # Only plot strains that have appeared  ################### THIS IS WRONG; MAKE IT "AT LEAST ONCE" RATHER THAN "SUM OF"
    tally = 0    
    for prob in p:
        tally += prob
        
    if (tally > 0):
        #plt.plot( w, p, label=lbl, linewidth='3' )
        plt.scatter( w, p, label=lbl, color=next(colors), s=400 )


plt.legend(bbox_to_anchor=(1.03, 1), loc=2, borderaxespad=0.)

plt.show()
#figureName = "MPMP_cutoff_"+str(CUTOFF)+".png"
#plt.savefig( join(path,figureName) )
#plt.savefig( join(path,figureName), bbox_inches='tight')

################################################################################################




exit(0)

