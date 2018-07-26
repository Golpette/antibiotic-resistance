#
# Script to produce fitness landscape in a given environment (i.e. [ciprofloxacin])
#
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

# Hashmap of genotypes and their fitnesses (wrt wildtype).
FITNESS = { "00000":1.0,
       "10000":1.01, "01000":0.99, "00100":0.99, "00010":0.83, "00001":0.91,      
       "11000":0.97, "10100":0.98, "10010":0.86, "10001":0.95, "01100":1.02, "01010":0.83, "01001":0.92, "00110":0.83, "00101":0.91, "00011":0.82,
       "11100":1.01, "11010":0.79, "11001":0.93, "10110":0.86, "10101":0.92, "10011":0.78, "01110":0.84, "01101":0.95, "01011":0.60, "00111":0.82,
       "11110":0.89, "11101":0.94, "11011":0.66, "10111":0.71, "01111":0.72,
       "11111":0.68
       }




# Read in antibiotic concentration.dat
path = join(".",str(sys.argv[1]))
antibifile = path+"/antibiotic.dat"



f = open(antibifile,'r')



for line in f:

    columns = [ float(x) for x in line.split() ]
    CONC =  columns[1]
    # string rep for file names
    cipcon = str( round(CONC,1) )



    ### Now make both types of plot for each well
    TRUE_FITNESS = {}
    for gen in all_genotypes:
        genome_label = "".join( ["%d" %x for x in gen] )    
        death = ( CONC/MICS[genome_label] )**2
        true_fit=0.0
        if ( death >= 1.0 ):
            true_fit = 0
        else:
            true_fit = FITNESS[genome_label] * ( 1 - death )
        TRUE_FITNESS[genome_label] = true_fit

    ##print TRUE_FITNESS



    # Get most fit genotype in this environment
    bestFitness = 0
    bestGeno = ""
    for key in TRUE_FITNESS:
        if TRUE_FITNESS[key] > bestFitness:
            bestFitness = TRUE_FITNESS[key]
            bestGeno = key

    ##print ""
    ##print bestGeno, bestFitness





    ####### Make plot of "fitness landscape" histogram #######
    genotypeInt = []
    genotypeStrings = []
    g_int=0
    for gen in all_genotypes:
        g_int +=1 
        genotypeInt.append( g_int )
        genome_label = "".join( ["%d" %x for x in gen] )   
        genotypeStrings.append( genome_label )
    trueFitnesses = []
    for g in genotypeStrings:
        trueFitnesses.append( TRUE_FITNESS[g] )

    plt.clf()
    title = "Ciprofloxacin = "+cipcon+" ng/ml"
    plt.title(title, fontsize=60)
    plt.bar( genotypeInt, trueFitnesses, align='center', color='olivedrab')
    plt.xticks(rotation=70, fontsize=20)
    plt.xticks( genotypeInt, genotypeStrings)
    plt.yticks( fontsize=20 )
    plt.ylabel( "Fitness", fontsize=40 )
    plt.ylim((0,1.2))
    text = "g* = "+bestGeno
    plt.text(15, 1.1, text, fontsize=32)
    #plt.show()
    figure = plt.gcf() # get current figure
    figure.set_size_inches(16, 12)
    # when saving, specify the DPI
    filename2 = path+"/histo_"+cipcon+".png"
    plt.savefig(filename2,dpi=100)








    plt.clf()

    #######  Make plot like Claudia Bank's "fitness landscape"  #######
    # lists of x,y coords and labels
    xp = []
    yp = []
    lab = []

    # dictionary of coords of all genotypes 
    coords = {}

    for gen in all_genotypes:
        # count # drivers
        cdriv = 0
        for g in gen:
            if( g==1 ):
                cdriv += 1
    
        genome_label = "".join( ["%d" %x for x in gen] )   

        xp.append( cdriv )
        yp.append( TRUE_FITNESS[ genome_label] )
        if(TRUE_FITNESS[ genome_label]>0):
            lab.append( genome_label )
        else:
            lab.append( "" )

        # store coords of each point so we can draw lines later
        coords[genome_label]=(cdriv,TRUE_FITNESS[genome_label])

    plt.scatter(xp, yp)

    for i, txt in enumerate( lab ):
        plt.annotate(txt, (xp[i],yp[i]))

    # connect all mutations
    # draw vertical line from (70,100) to (70, 250)
    #plt.plot([70, 70], [100, 250], 'k-', lw=2)

    # Enumerate all reactions with DFS
    rxn_list = []
    no_drivers = 0
    rxn_list = enum_rxns( [0,0,0,0,0], rxn_list, no_drivers )
    # Print line between for all possible reactions
    for (r1, r2) in rxn_list:
        rxn1 = list(r1)
        rxn2 = list(r2)
        g1 = "";  g2 = "";
        for pm in r1:
            g1 += str(pm)
        for pm in r2:
            g2 += str(pm)
    
        x1,y1 = coords[g1]
        x2,y2 = coords[g2]
        plt.plot([x1, x2], [y1, y2], 'k-', lw=0.2)

    plt.xticks(fontsize=20)
    plt.yticks( fontsize=20 )
    plt.ylabel( "Fitness", fontsize=40 )
    plt.xlabel( "Number of mutations", fontsize=40 )
    plt.ylim((0,1.1))
    #plt.ylim(ymax=1.1)
    #plt.show()
    filename3 = path+"/lines_"+cipcon+".png"
    plt.savefig(filename3)













