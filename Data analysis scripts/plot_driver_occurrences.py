#################################################################################################
#
# Program to plot in which well each driver mutation becomes measurable (abundance<CUTOFF)
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


# Input directory containg the "well_x_...dat" files
path = join(".", str(sys.argv[1]))

# Define minimum abundance of PM to be experimentally detectable (about 10% is likely)
CUTOFF = 0.1

# Define genomic position of driver PMs
DRIVERS = [0,1,2,3,4]

# List of wells we will get stats for
wells = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]



num_experiments = 5    # Set by hand here so that we don't always print hundreds of images



for exp in range(0,num_experiments):

    # Store (well, driver PM) tuples for each experiment
    well_driver = []

    for w in wells:

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
                
                # Read in line as list of numbers (rank, PM_location, abundance) and store only
                # driver PMs with abundance > CUTOFF
                columns = [ float(x) for x in line.split() ] 
                if ( columns[2] > CUTOFF and columns[1] in DRIVERS ):
                    well_driver.append(  (w, columns[1])   ) 



    ###################### PLOTTING DATA FOR EACH EXPERIMENT #######################################
    fig = plt.figure()
    fig.set_size_inches(14, 8)
    title = "Driver PMs present in wells (cutoff = "+str(CUTOFF)+")"
    fig.suptitle(title, fontsize=32)
    plt.xlabel('Well', fontsize=26);  
    #plt.ylabel('Driver', fontsize=22)
    plt.xlim((0,25));  plt.ylim((-1,5))
    plt.tick_params(axis='both', which='major', labelsize=24)

    axes = fig.add_subplot(111)
    a=axes.get_yticks().tolist()
    a[0]=''; a[1]='gyrA1';  a[2]='gyrA2';  a[3]='parC';  a[4]='marR';  a[5]='acrR';  a[6]='';
    axes.set_yticklabels(a)

    # Make column scatter plot of "well vs #PMs" for all detectable PMs
    x=[];  y=[];
    for (a,b) in well_driver:
      x.append(a);   y.append(b);
    plt.scatter(x, y, color='green', marker='s', s=800)

    #plt.show()
    figureName = "Experiment_"+str(exp)+"_cutoff_"+str(CUTOFF)+".png"

    # Only save figure if experiment exists
    if len(required_files) > 0:
        plt.savefig( join(path,figureName) )
        #plt.savefig( join(path,figureName), bbox_inches='tight')
    ################################################################################################







exit(0)

