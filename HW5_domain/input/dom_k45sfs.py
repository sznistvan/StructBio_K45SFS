#!/usr/bin/env python
#This notebook is a guide on how to start a Python3 code that implements the pseudodode provided for the domain assignment task. 
#The test file is "1GPZ.pdb", which should be in the same folder as this notebook file (domain.ipynb). 
#I will provide detailed help on how to read the file in but will not show an entire working code in this notebook - it is your task to put the pieces together and complete the code. 
#Note that he pieces of code shown might (should :-) appear in the full working program in an order different from that shown here. 
#You are of course encouraged to prepare your own code in a style you prefer.

#During PDB file processing, x,y and z coordinates of each CA atom will be stored in arrays. The array indices, ranging from 0 to the number of residues found, are not expected to correspond to the residue numbering in the PDB file (in our case, starting from 290), thus, the residue numbers as shown in the PDB file will also be stored. The data structures will be as follows:


import numpy
import matplotlib.pyplot as mt
import sys

# storing coordinates in arrays ranging from 0 to
# the number of residues
coordx = []
coordy = []
coordz = []
# residue numbers will also be stored as numbering in PDB files
# rarely starts with 1, residues can be missing etc.
pdb_resnum = []


# We will also have to initialize two other variables, the number of residues that will be increased after reading and processing each relevant line (containing coordinates of CA atoms of residues in chain A), and a cutoff value to define contacts. The cutoff value in the example will be set to 8, meaning that two residues will be considered to be in contact with each other if the distance between their CA atoms is no more than 8 Angstroms.



# cutoff value for CA-CA distances
cutoff = 8
# number of residues that have been read in
nres = 0

protein_id =""
error = False


# Below is an example of a loop reading the coordinates of CA atoms in chain A from the PDB file. Note that the PDB file reading shown below is a very simple one. Each field will be read based on their positions in the line (that follows from the specification of the line), but only the relevant fileds will be read and from 'ATOM' lines.
# For the test file, only chain 'A' will be processed.



# defining the name of the file. Should be in the same folder as the program/notebook.
if(len(sys.argv)>=2):
    try:
        pdbfile = open(sys.argv[1],"r")
        if(pdbfile != None):
            print("ok")
            for line in pdbfile:
                if line.startswith("HEADER"):
                    protein_id=str(line[62:66])
                # only ATOM lines will be considered
                if line.startswith("ATOM"):
                    # obtaining the atom name and getting rid of spaces
                    atom=line[12:15].replace(" ","")
                    # the chain identifier is 1 character
                    chain=line[21:22]
                    # processing only lines of CA atoms in chain E
                    if atom == "CA" and chain == sys.argv[2]:
                        # the append function adds a value to the array. The substring at given positions of the line
                        # should be converted to a floating-point number.
                        coordx.append(float(line[30:37]))
                        coordy.append(float(line[38:45]))
                        coordz.append(float(line[46:53]))
                        # storing residue numbers, they are integers 
                        pdb_resnum.append(int(line[22:27]))
                        # we have appended a value to each array, they should have equal size and the data
                        # for the same residue should be accessible at the same index that corresponds to nres 
                        # (the first index is zero and we set nres to zero above). Now we increase nres to keep track.
                        nres=nres+1
    except FileNotFoundError:
        print("ERROR! PDB file doesn't exist!")
        error=True
else:
    print("ERROR! Usage: python {} PROT.pdb CHAIN".format(sys.argv[0])) 
    error = True


# Now we know how many residues we have, we can now initialize the contact matrix. It will have a size of nres x nres and we will fill it with zeroes, the code below does this.
# Initialization: filling the contact matrix with zeroes
contact_matrix = [[0 for h in range(nres)] for k in range(nres)]


# Now we have to fill the contact matrix by calculating the distance of every residue pair i,j. We only have to calculate the upper triangle part of the matrix where 0 < i < nres and 0 <= j < i. A contact will be added by setting the matrix value to 1 when the distance between the CA atoms is not larger than the cutoff. The distance is an Euclidean distance to be calculated using Pythagoras' formula.
# I recommend to make the code "foolproof" by filling the matrix in a symmetric fashion, i.e. ensuring that



for i in range(1, nres):
    for j in range(0, i): 
        distance=((coordx[i]-coordx[j])*(coordx[i]-coordx[j]))+((coordy[i]-coordy[j])*(coordy[i]-coordy[j]))+((coordz[i]-coordz[j])*(coordz[i]-coordz[j]))
        distance=numpy.sqrt(distance)

        if distance <= cutoff:
            contact_matrix[j][i]=1
            contact_matrix[i][j]=contact_matrix[j][i]
        #print("i",i," j",j,"   distance",distance," matrix",contact_matrix[i][j])
            
# As we already have the full matrix this will not use extra storage and will make the code work even if we forget whether we calculated the upper or he lower triangle in the first place :-)
# Now we have all the contacts. We have now to iterate over all possible chain cuts and calculate the number of intra- and interdomain contacts. For this, the definition of functions is recommended:


def contacts_intradomain(p, q):
    contacts_i=0
    for r in range(p,q):
        for s in range(p,r):
            contacts_i+=contact_matrix[r][s]
            
    return contacts_i
    
def contacts_interdomain(p1, p2, q1, q2):
    contacts_e=0
    for r in range(p1,q1):
        for s in range(p2,q2):
            contacts_e+=contact_matrix[r][s]

    return contacts_e

to_plot=[]
writefile = open("output.txt","a")
for c in range (1, nres-1):
    intradomain_A=contacts_intradomain(0,c)
    intradomain_B=contacts_intradomain(c+1,nres)
    interdomain_AB=contacts_interdomain(0,c+1,c,nres)
    d=intradomain_A*intradomain_B/(interdomain_AB*interdomain_AB)
    print("A:", intradomain_A," B:",intradomain_B,"  AB:",interdomain_AB, "{} \t {}".format(pdb_resnum[c],d))
    to_plot.append(d)
    outp = "A:"+ str(intradomain_A)+" B:"+str(intradomain_B)+"  AB:"+str(interdomain_AB)+ "{} \t {}".format(pdb_resnum[c],d)+"\n"
    writefile.write(outp)


mt.plot(to_plot,color="red")
mt.title("Domain identification ({})".format(protein_id))
mt.xlabel("Residues")
mt.ylabel("IntraA * IntraB / InterAB^2")
if(not error):
    mt.show()




