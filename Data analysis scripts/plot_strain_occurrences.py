#################################################################################################
#
# Program to plot abundances of different mutant strains across wells at end of experiments
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


#################################################################################################
num_experiments = 10   # no. of experiments to plot

num_strains = 32


# Input directory containg the "strain__x_...dat" files
path = join(".", str(sys.argv[1]))

# Define minimum abundance of strain to be experimentally detectable 
#   (doesn't really make sense as before when we were talking about detecting specific mutated genes)
CUTOFF = 0.1

# List of wells we will get stats for
wells = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]
#################################################################################################



## QUICK FIX: I was naming strains with integers before. Stupid. Want to print actual genotype as label
strainMap = {0:'00000', 1:'10000',2:'01000',3:'00100',4:'00010',5:'00001',
             6:'11000',7:'10100',8:'10010',9:'10001',10:'01100',11:'01010',12:'01001',13:'00011',
             14:'11100',15:'11010',16:'11001',17:'10110',18:'10101',19:'01110',20:'01101',21:'10011',22:'01011',
             23:'11011',24:'11110',25:'11101',26:'10111',27:'01111',28:'11111',
             29:'00110',30:'00101',31:'00111'}
#### This is hardcoded in the main program and is the order they appear in Marcusson2009 (+3 extras at end)




for exp in range(0,num_experiments):

    # Store (strain, well, number) tuples for each experiment
    experiment_data = []

    for w in wells:

        # Root of files relating to current well
        search = "strains_w"+str(w)+"___(.*)_"+str(exp)+".dat"

        # Get file for given experiment and well 
        #required_files = [f for f in listdir('.') if( isfile(f) and re.search(search, f) ) ]
        required_files = [f for f in listdir(path) if( isfile(join(path,f)) and re.search(search, f) ) ] 

        if len(required_files) > 1: 
            print "ERROR: more than 1 file found for given experiment and well."
            print required_files, "###"
        
        for filename in required_files:
            
            f = open( join(path,filename), 'r' )
            for line in f:
                
                # Read in line as list of numbers (strain, number) and store only
                # driver PMs with abundance > CUTOFF
                columns = [ float(x) for x in line.split() ] 
                if ( columns[1] > 0 ):
                    # Store tuples of (strain, well, number)
                    experiment_data.append(  (columns[0], w, columns[1])   )
                    ##well_driver.append(  (w, columns[1])   ) 



    ###################### PLOTTING DATA FOR EACH EXPERIMENT #######################################
    fig = plt.figure()
    fig.set_size_inches(14, 8)
    title = "Final strain abundances"
    fig.suptitle(title, fontsize=32)
    plt.xlabel('Well', fontsize=26);  
    plt.ylabel('Population / K', fontsize=26)
    plt.xlim((0,25));  plt.ylim((0,1.1))
    plt.tick_params(axis='both', which='major', labelsize=22)

    #axes = fig.add_subplot(111)
    #a=axes.get_yticks().tolist()
    #a[0]=''; a[1]='gyrA1';  a[2]='gyrA2';  a[3]='parC';  a[4]='marR';  a[5]='acrR';  a[6]='';
    #axes.set_yticklabels(a)

    for s in range(num_strains):
        # Make column scatter plot of "well vs #PMs" for all detectable PMs
        #strain=[]; 
        w=[];  N=[];
        for (a,b,c) in experiment_data:
            if( a == s ):
                w.append(b);   N.append(c);
                #plt.scatter(w, N)
                # Only plot if at least 1 well is non-zero
        tot_n = 0
        for n in N:
            tot_n += n 
        if( tot_n > CUTOFF ):
            #lbl = "Strain "+str(s)
            lbl = strainMap[s]
            plt.plot(w,N,linewidth='3',label=lbl)

    #plt.legend(bbox_to_anchor=(0.95, 1), loc=2, borderaxespad=0.)

    #plt.show()
    #figureName = "Experiment_"+str(exp)+"_cutoff_"+str(CUTOFF)+".png"
    #figureName = "Strains_Experiment_"+str(exp)+".png"

    #plt.savefig( join(path,figureName) )

    # Only save figure if experiment exists
    if len(required_files) > 0:

        plt.legend(bbox_to_anchor=(0.95, 1), loc=2, borderaxespad=0.)
        figureName = "Strains_Experiment_"+str(exp)+".png"

        plt.savefig( join(path,figureName) )
    #    #plt.savefig( join(path,figureName), bbox_inches='tight')
    ################################################################################################







exit(0)

