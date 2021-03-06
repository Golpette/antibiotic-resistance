#
# Script to take ordered list of driver mutations ("sequence") and convert
# it to a mutational pathway (of form: 00000 00010 00011 etc).
# Also produces weighted graphviz file.
#
import os
from os import listdir
from os.path import isfile, join
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





# Edge widths are frequency * this arbitrary scaling
ARBITRARY_SCALING = 100.0
DRIVERS = 5


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






path = join(".", str(sys.argv[1]))
f = open( path,'r')


fileout = "MPMPs.dat"
fout = open(fileout, 'w')

### Read in sequences (i.e. final genotypes) "4 3 2 1" from Bartek's program
###  and print mutational pathway: "00000 00001 00011..." etc
for line in f:

    genome = ["0" for i in range(DRIVERS)]

    # Read in line: "number" and trajectory ("4 3 2 1" etc.)
    columns = [ int(x) for x in line.split() ]

    # Getand print number that followed this trajectory
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

   
#exit(0)



########## Store each mutation and occurrence in dictionary ################

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


print mutation_counts
print "counter = ", counter

fin.close()
###########################################################################


#exit(0)


############  Make graphviz file  #########################################

# Enumerate all reactions
rxn_list = []
no_drivers = 0
rxn_list = enum_rxns( [0,0,0,0,0], rxn_list, no_drivers )


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

        # Make string from integer list
        genome_label = "".join( ["%d" %x for x in g] )

        # Get color from MICS hashmap !!!NOTE: Each color scheme has a different range; e.g. "orrd9" is [1,9])!!!
        nodecolor = int( norm( MICS[genome_label] )*8 )+1

        # Write node object: label, height, position, width
        gv.write( genome_label + " [height=0.3,pos=\""+str(Xpos)+","+str(Ypos)+"!\", width=1.0, style=filled, fillcolor="+str(nodecolor)+"];\n"   )



# Print all reactions with no edges
for (r1, r2) in rxn_list:
    rxn1 = list(r1)
    rxn2 = list(r2)
    g1 = "";  g2 = "";
    for pm in r1:
        g1 += str(pm)
    for pm in r2:
        g2 += str(pm)
    gv.write( g1 + " -- " + g2 + " [penwidth=0];\n" )


# Print modified reactions with weighted edges
for key in mutation_counts:
    fraction = float(mutation_counts[key])/float(counter)
    # scale line width ARBITRARILY with fraction; will depend on number of mutations (i.e. counter)
    line_width_zz = fraction*ARBITRARY_SCALING
    #format to avoid scientific notation:
    line_width = '{:.9f}'.format( line_width_zz )    
    gv.write( key+" [penwidth="+str(line_width)+"];\n" )
    #gv.write( key+" [label=\""+str(fraction)+"\", penwidth="+str(line_width)+"];\n" )


gv.write( "}" )
gv.close()


# Make image
os.system("dot -s -Kneato -n -Tpng -o OUTPUT.png graphviz_file.dat")

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


