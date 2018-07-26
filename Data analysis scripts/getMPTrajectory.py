
########## Get most probable trajectory observed from MPMP.dat ############
fin_2 = open("MPMPs.dat", 'r')

trajectory_counts = {}

for line in fin_2:

    # Read in columns as strings
    columns = [ s for s in line.split() ]
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







while( trajectory_counts ):
    highest = 0
    highestkey = ""
    for key in trajectory_counts:
        if trajectory_counts[key] > highest:
            highest = trajectory_counts[key]
            highestkey = key
    print highestkey + " " + str(highest)
    del trajectory_counts[highestkey]



exit(0)





'''
# Get most observed trajectory
highest = 0
highestkey = ""
for key in trajectory_counts:
    if trajectory_counts[key] > highest:
        highest = trajectory_counts[key]
        highestkey = key


#Split this trajectory into it's reactions:
genos = highestkey.split()

rxns_in_MPTraj = []
for i in range( len(genos)-1 ):
    rxn = genos[i]+" -- "+genos[i+1]
    rxns_in_MPTraj.append(rxn)

print genos
print rxns_in_MPTraj

exit(0)
'''




###########################################################################

