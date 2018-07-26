import os
import sys
from os.path import join


# Read in directory where "sequence_)wellX...dat" files are
#dirpath = join(".", str(sys.argv[1]) )
dirpath = str(sys.argv[1]) 


for i in range(24):

    command = "python construct_geno_from_seq_MPEP_changeSeqFormat_MPTraj.py " + dirpath + " " + str(i)

    #print command
    #exit(0)

    os.system( command )

