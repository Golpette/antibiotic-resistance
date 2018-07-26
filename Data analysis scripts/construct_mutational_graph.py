#
# Script to take ordered list of driver mutations ("sequence") and convert
# it to a mutational pathway (of form: 00000 00010 00011 etc).
# Also produces weighted graphviz file.
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



def enum_rxns( geno, rxn_list, no_drivers ):
    ''' DFS algorithm to enumerate all reactions '''
    # if we give [0,0,0,0,0] as start genotype

    for x in range( len(geno) ):
        newlist = list(geno)
        newlist[ x ] = 1
        tot=0
        for n in newlist:
            tot += n
            if( tot > no_drivers ):
                # Only want individual reactions, not every single 
                # mutational pathway
                if( (geno,newlist) not in rxn_list ):
                    rxn_list.append( (geno,newlist) )
                    no_drivers += 1
                    enum_rxns( newlist, rxn_list, no_drivers )
                    no_drivers -= 1
    return rxn_list






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




# Edge widths are frequency * this arbitrary figure
ARBITRARY_SCALING = 100.0
DRIVERS = 5


# Input directory containing the "strain_wx___...dat" files
dirpath = join(".", str(sys.argv[1]))
wellsearch = str(sys.argv[2])   # i.e. "23"  (or "*" or "!23" if possible?)



# Cat all desired wells into single file: 
#   (need to escape brackets here since they're going into terminal command)
collatedfile_zz = join(dirpath,"_allSequences\(well"+wellsearch+"\).dat")
command = "cat "+dirpath+"sequences_\(well"+wellsearch+"\)_*.dat > "+collatedfile_zz 
os.system(command)
# Actual collated filename as string varibale (i.e. no escaped characters)
collatedfile = join(dirpath,"_allSequences(well"+wellsearch+").dat")


# Root of files relating to current well
search = "sequences_\(well"+wellsearch+"\).*.dat"
# Make list of desired files 
required_files = [f for f in listdir(dirpath) if( isfile(join(dirpath,f)) and re.search(search, f) ) ]




########## Count number of experiments (by counting well0 files...) ##############
# Root of files relating to first well
search_well0 = "sequences_\(well0\)"+".*.dat"
# Make list of desired files 
well0_files = [f for f in listdir(dirpath) if( isfile(join(dirpath,f)) and re.search(search_well0, f) ) ]

NUM_EXPERIMENTS = len( well0_files )
#################################################################################




####### Check all files to calculate prob: P(genotype present in well > cutoff) ####################

# store which genotypes appear and in how many experiments
geno_occurrences = {}

for filename in required_files:
    
    ff = open( join(dirpath,filename), 'r')

    this_experiment = {}
    this_population = 0
    
    for line in ff:

        #(first entry is the number of cells)
        mutatedsites = [ int(x) for x in line.split() ]
        numberofcells = mutatedsites.pop(0)
        this_population += numberofcells

        # create final genotype 
        g = ""
        final_g = ["0" for i in range(DRIVERS)]  
        for driver in mutatedsites:
            final_g[ driver ] = "1"
            g = "".join(final_g)
        if g == "":
            # i.e. no mutated sites
            g = "00000"
        if g in this_experiment:
            this_experiment[g] += numberofcells
        else:
            this_experiment[g] = numberofcells 

    # count genotypes present above 0.1 abundance
    cutoff = float(this_population/10.0)

    for gen in this_experiment:
        if this_experiment[gen] > cutoff:
            if gen in geno_occurrences:
                geno_occurrences[gen] += 1
            else:
                geno_occurrences[gen] = 1

#print all
#print geno_occurrences
#exit(0)

#########################################################################




### Read in sequences (i.e. final genotypes) "4 3 2 1" from Bartek's program
###  and print mutational pathway: "00000 00001 00011..." etc

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




########## Store each mutation and occurrence in dictionary ################
# Read in the MPMP file
fin = open( fileout, 'r' )

mutation_counts = {} 

counter = 0  # used to weight graph edges (counts *absolute* number of mutations)

for line in fin:
    
    # Read in columns as strings
    columns = [ s for s in line.split() ]

    # Number that followed this trajectory
    number = int(columns.pop(0))

    # Add all mutations and counts to dictionary
    c = 0
    for s in columns:
        if( c < len(columns)-1 ):
            mut_string = s+" -- "+columns[c+1]
            c += 1
            counter += number
            if mut_string in mutation_counts:
                mutation_counts[ mut_string ] += number
            else:
                mutation_counts[ mut_string ] = number

#print mutation_counts
#print "counter = ", counter

fin.close()
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





############  Make graphviz file  #########################################

# Logarthmic normalization for MIC-colouring of nodes
norm = matplotlib.colors.LogNorm(vmin=32.0, vmax=64000.0) 


gv = open( "graphviz_file.dat", 'w' )

# Print preamble
gv.write( "graph {\n" )
gv.write( "graph [bb=\"0,0,1200,500\"];\n" )
gv.write( "node [label=\"\N\", fontname=\"times bold\", fontsize=16, penwidth=0, colorscheme=\"orrd9\"];\n" )



# Print all nodes and their *fixed* positions
Xmax = 1000;  Ymax = 600;  Yinc = float(Ymax)/11.0;
five_choose_n = [1,5,10,10,5,1]

for nds in range(DRIVERS+1):
    temp_list = []
    for gen in all_genotypes:

        totnds=0  # count drivers
        for n in gen:
            totnds += n

        if( totnds == nds ):
            temp_list.append( gen )
    cc=0        
    for g in temp_list:
        cc += 1
        # divide x-range into no. of drivers+1 (start from 0)
        Xint = float(Xmax)/float( DRIVERS + 2.0 )
        Xpos = (1+nds)*Xint
        Ypos = cc*Yinc

        # Make string from intger list
        genome_label = "".join( ["%d" %x for x in g] )

        # Get color from MICS hashmap !!!(Each color scheme has different ranges; "orrd9" is [1,9])!!!
        nodecolor = int( norm( MICS[genome_label] )*8 )+1

        # Write node object: label, height, position, width
        gv.write( genome_label + " [height=0.3,pos=\""+str(Xpos)+","+str(Ypos)+"!\", width=1.0, style=filled, fillcolor="+str(nodecolor)+"];\n"   )
        #gv.write( genome_label + " [height=0.1,pos=\""+str(Xpos)+","+str(Ypos)+"!\", width=0.1];\n"   )


# Enumerate all reactions with DFS
rxn_list = []
no_drivers = 0
rxn_list = enum_rxns( [0,0,0,0,0], rxn_list, no_drivers )
# Print all reactions with no edges (seems this makes edges clearer)
for (r1, r2) in rxn_list:
    rxn1 = list(r1)
    rxn2 = list(r2)
    g1 = "";  g2 = "";
    for pm in r1:
        g1 += str(pm)
    for pm in r2:
        g2 += str(pm)
    #gv.write( g1 + " -- " + g2 + ";\n" )
    gv.write( g1 + " -- " + g2 + " [penwidth=0];\n" )


# Print modified reactions with weighted edges
for key in mutation_counts:
    fraction = float(mutation_counts[key])/float(counter)
    line_width_zz = fraction*ARBITRARY_SCALING
    #format to avoid scientific notation:
    line_width = '{:.9f}'.format( line_width_zz )   

    # Write reaction in different colour if it belongs in MPTraj
    if (key in rxns_in_MPTraj ): 
        gv.write( key+" [color=\"forestgreen\", penwidth="+str(line_width)+"];\n" )
    else:
        gv.write( key+" [penwidth="+str(line_width)+"];\n" )


   #gv.write( key+" [label=\""+str(fraction)+"\", penwidth="+str(line_width)+"];\n" )



# Print fake nodes as legend, ordered with most likely endpoint
#sumcount=0
#for key in geno_occurrences:
#    sumcount += geno_occurrences[key]
#
# Probability genotype is present in final well > cutoff
geno_prob = {}
for key in geno_occurrences:
    #pr = 100*float(geno_occurrences[key])/float(sumcount)
    pr = 100*float(geno_occurrences[key])/float(NUM_EXPERIMENTS)
    geno_prob[key] = round(pr,1)
# Print highest first
gv.write( "9999 [label=\"P(genotype abundance > 10%)\", fontname=\"helvetica\", fontsize=18, pos=\"0,"+str(Ymax-75)+"!\"];\n")
#Controls arbitrary numeric label and Y-position of fake node
nodelabel=30
# Keep printing and removing most likely geno until none left
while( geno_prob ):
    highest = 0
    highestkey = ""
    for key in geno_prob:
        if geno_prob[key] > highest:
            highest = geno_prob[key]
            highestkey = key
    gv.write( str(nodelabel)+" [label = \""+highestkey+" = "+str(geno_prob[highestkey])+"%\", fontname=\"helvetica\", fontsize=18, pos=\"0,"+str(Ymax-75-nodelabel)+"!\"];\n"  )
    nodelabel += 30
    del geno_prob[highestkey]




gv.write( "}" )
gv.close()


# Make image
os.system("dot -s -Kneato -n -Tpng -o OUTPUT_"+wellsearch+".png graphviz_file.dat")

##########################################################################













"""
def enum_genotypes( geno, geno_list, no_drivers ):
    ''' DFS algorithm to enumerate all genotypes '''
    # if we give [0,0,0,0,0] as start genotype

    for x in range( len(geno) ):
        newlist = list(geno)
        newlist[ x ] = 1
        tot=0
        for n in newlist:
            tot += n
        if( tot > no_drivers ):
            if( newlist not in geno_list ):
                geno_list.append( newlist )
                no_drivers += 1
                enum_genotypes( newlist, geno_list, no_drivers )
                no_drivers -= 1
    return geno_list
"""


