#
# Script to plot rank-ordered cumulative probability of all trajectories
#   found in all cells in a given well, with specified end-point genotype
#
import os
from os import listdir
from os.path import isfile, join
# regexp; need search() method for substring in filename
import re
# arguments from command line
import sys
# for logarithmic normalizer color map
import matplotlib
import pylab as plt




# This is the order they will be placed when creating the graphviz file (is there a more optimal layout?)
all_genotypes = [ [0,0,0,0,0],
                [1,0,0,0,0], [0,1,0,0,0], [0,0,1,0,0], [0,0,0,1,0], [0,0,0,0,1],      
                [1,1,0,0,0], [1,0,1,0,0], [1,0,0,1,0], [1,0,0,0,1], [0,1,1,0,0], [0,1,0,1,0], [0,1,0,0,1], [0,0,1,1,0], [0,0,1,0,1], [0,0,0,1,1],
                [1,1,1,0,0], [1,1,0,1,0], [1,1,0,0,1], [1,0,1,1,0], [1,0,1,0,1], [1,0,0,1,1], [0,1,1,1,0], [0,1,1,0,1], [0,1,0,1,1], [0,0,1,1,1],
                [1,1,1,1,0], [1,1,1,0,1], [1,1,0,1,1], [1,0,1,1,1], [0,1,1,1,1],
                [1,1,1,1,1]
                ]
# Hashmap of genotypes and their MICs (these are 2 times the Marcusson2009 value)
MICS = { "00000":32,
       "10000":760, "01000":500, "00100":32, "00010":64, "00001":94,      
       "11000":760, "10100":2000, "10010":2000, "10001":1000, "01100":760, "01010":2000, "01001":760, "00110":64, "00101":94, "00011":250,
       "11100":64000, "11010":1500, "11001":1500, "10110":12000, "10101":6000, "10011":3000, "01110":1500, "01101":1000, "01011":3000, "00111":250,
       "11110":64000, "11101":64000, "11011":4000, "10111":16000, "01111":4000,
       "11111":64000
       }


DRIVERS = 5




# Input directory containing the "strain_wx___...dat" files
dirpath = join(".", str(sys.argv[1]))
wellsearch = str(sys.argv[2])   # i.e. "23"  (or "*" or "!23" if possible?)
desired_end_gen = str( sys.argv[3] )


# Cat all desired wells into single file: 
#   (need to escape brackets here since they're going into terminal command)
collatedfile_zz = join(dirpath,"_allSequences\(well"+wellsearch+"\).dat")
command = "cat "+dirpath+"sequences_\(well"+wellsearch+"\)_*.dat > "+collatedfile_zz 
os.system(command)
# Actual collated filename as string variable (i.e. no escaped characters)
collatedfile = join(dirpath,"_allSequences(well"+wellsearch+").dat")


# Root of files relating to current well
search = "sequences_\(well"+wellsearch+"\).*.dat"
# Make list of desired files 
required_files = [f for f in listdir(dirpath) if( isfile(join(dirpath,f)) and re.search(search, f) ) ]


NUM_EXPERIMENTS = len( required_files )








### Read in sequences (i.e. final genotypes) "4 2 3 1" from Bartek's program
###  and print mutational pathway: "00000 00001 00101..." etc

#path = join(".", str(sys.argv[3]))
f = open( collatedfile,'r')

fileout = "MPMPs.dat"
fout = open(fileout, 'w')

for line in f:

    genome = ["0" for i in range(DRIVERS)]

    # Read in line: "number" and trajectory ("4 2 1 3" etc.)
    columns = [ int(x) for x in line.split() ]

    # Get and print number that followed this trajectory
    number = columns.pop(0)
    fout.write( str(number) + " " )
    fout.write( "".join(genome) )  ############## i.e. initial "00000" #CAREFUL! MAYBE WE DON'T ALWAYS WANT TO START FROM WT ???
    fout.write(" ")

    for driver in columns:

        genome[ driver ] = "1"
        fout.write( "".join(genome) )
        fout.write(" ")

    fout.write("\n")

f.close()
fout.close()
##########################################################################






########## Make plot of cumulative probabilities of trajectories in MPMP.dat ############
fin_2 = open(fileout, 'r')

trajectory_counts = {}

for line in fin_2:

    # Read in columns as strings
    columns = [ s for s in line.split() ]

    # ONLY CONSIDER IF DESIRED_FINAL_GEN IS REACHED
    # Want only "trajectories" so ignore WT genotypes (i.e. entries with only "number 00000":
    cl = len(columns)
    if cl > 2 and columns[ cl-1 ]==desired_end_gen :

        # Number that followed this trajectory
        number = int(columns.pop(0))

        # Put remaining columns back together as trajectory
        traj = ""
        for s in columns:
            traj = traj + " " + s
 
        if traj in trajectory_counts:
            trajectory_counts[traj] += number
        else:
            trajectory_counts[traj] = number


# Sum total number of trajectories (i.e. cells across all experiments)
totalCells = 0
for key in trajectory_counts:
    totalCells += trajectory_counts[key]


# Get rank and probability of all trajectories
ranks = []
probs = []
pths = []

rnk = 0
prb = 0
while( trajectory_counts ):
    highest = 0
    highestkey = ""
    for key in trajectory_counts:
        if trajectory_counts[key] > highest:
            highest = trajectory_counts[key]
            highestkey = key

    rnk += 1
    ranks.append( rnk )
    prb += float( trajectory_counts[highestkey] ) / float( totalCells )
    probs.append( prb )
    ##print highestkey + " " + str(highest)

    # Store each trajectory as a string 
    traj=""
    # Split this trajectory into it's reactions:
    genos = highestkey.split()
    for i in range( len(genos)-1 ):
       #traj += str(genos[i])+" -- "+str(genos[i+1])
        traj += genos[i]+"--"
    traj += genos[ len(genos)-1 ]
    pths.append( traj )

    del trajectory_counts[highestkey]


# Make plot
plt.xlabel("Rank", fontsize=26)
plt.ylabel("Cumulative probability", fontsize=26)
#plt.title()
plt.xticks(fontsize=20);    plt.yticks(fontsize=20)
plt.xlim((0,len(ranks)+1))
plt.ylim((0,1.05))

plt.plot( ranks, probs, linewidth=2, marker='o', markersize=10 )

#Put top 5 trajectories on plot
plt.text(0.4*len(ranks),0.6,"1: "+pths[0],fontsize=22)
plt.text(0.4*len(ranks),0.55,"2: "+pths[1],fontsize=22)
plt.text(0.4*len(ranks),0.5,"3: "+pths[2],fontsize=22)
plt.text(0.4*len(ranks),0.45,"4: "+pths[3],fontsize=22)
plt.text(0.4*len(ranks),0.4,"5: "+pths[4],fontsize=22)

plt.show()



###########################################################################





















































########## Get most probable TRAJECTORY observed from MPMP.dat ############
# Store the rxns present in the most common trajectory
rxns_in_MPTraj = []

fin_2 = open(fileout, 'r')

trajectory_counts = {}

for line in fin_2:

    # Read in columns as strings
    columns = [ s for s in line.split() ]

    # Want only "trajectories" so ignore WT genotypes
    #  (i.e. entries with only "number 00000":
    if len(columns) > 2:

        # Number that followed this trajectory
        number = int(columns.pop(0))

        # Put remaining columns back together as trajectory
        traj = ""
        for s in columns:
            traj = traj + " " + s
 
        if traj in trajectory_counts:
            trajectory_counts[traj] += number
        else:
            trajectory_counts[traj] = number



#while( trajectory_counts ):
#    highest = 0
#    highestkey = ""
#    for key in trajectory_counts:
#        if trajectory_counts[key] > highest:
#            highest = trajectory_counts[key]
#            highestkey = key
#    print highestkey + " " + str(highest)
#    del trajectory_counts[highestkey]


# Get "single" most observed trajectory
highest = 0
highestkey = ""
for key in trajectory_counts:
    if trajectory_counts[key] > highest:
        highest = trajectory_counts[key]
        highestkey = key
# Split this trajectory into it's reactions:
genos = highestkey.split()
# Store rxns present in most probably trajectory (MPTraj)
for i in range( len(genos)-1 ):
    rxn = genos[i]+" -- "+genos[i+1]
    rxns_in_MPTraj.append(rxn)



###########################################################################





