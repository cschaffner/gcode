"""
Computing the molecular structure matrix, an invention of Peter van der Gulik
It is a distance measure of the 20 amino acids. For two amino acids a1 and a2,
we use the SMSD (Small Molecule Subgraph Detector) solver from
http://www.ebi.ac.uk/thornton-srv/software/SMSD/
to find their maximum common subgraph (MCS), called "sub".

Then, we count the number of bonds that remain in a1 minus sub
and add them to the number of bonds remaining in a2 minus sub. 

Output is the lower triangular matrix with the pairwise distances, in latex-compatible format.

python script written by Christian Schaffner, c.schaffner@uva.nl
19 July 2012
"""

import os
import subprocess
import re
from numpy import arange

class Bond:
    """ Bonds have three properties: BeginAtomIdx, EndAtomIdx, BondType (1,2) """
    def __init__(self,i1,i2,typ):
        self.BeginAtomIdx=i1
        self.EndAtomIdx=i2
        self.BondType=typ        
        if self.EndAtomIdx == None or self.BondType<0 or self.BondType>2:
            raise

class Molecule:
    def __init__(self,molfile):
        """ imports bonds from a molfile """
        self.Bonds=[]        
        f = open(molfile, 'r')
        f.readline()
        f.readline()
        f.readline()
        line4=f.readline()
        # parse line 4 into elements
        el4=re.findall('\w+',line4)
        if el4[-1]!='V2000':
            print "mol version unknown: {0}".format(line4)
            raise
        nr_atoms=int(el4[0])
        nr_bonds=int(el4[1])
        if nr_atoms==0:
            print 'molecule has no atoms'
            raise
        if nr_bonds==0:
            print 'molecule has no bonds'
        # forward to bonds:
        for i in arange(nr_atoms):
            f.readline()
        # read in bonds
        for i in arange(nr_bonds):
            b_nrs=re.findall('\d+',f.readline())
            newBond=Bond(int(b_nrs[0]),int(b_nrs[1]),int(b_nrs[2]))
            self.Bonds.append(newBond)
            
    def GetBonds(self):
        for b in self.Bonds:
            yield b

def bonds_notin(m,match):
# returns the number of bonds that are in m, but not in the list match

#    print match    
    nr_bonds=0
    for bnd in m.GetBonds():
        if bnd.BeginAtomIdx in match and bnd.EndAtomIdx in match:
            continue
#        print 'bond nr {0}: {1} - {2}, {3}'.format(nr_bonds,bnd.GetBeginAtomIdx(),bnd.GetEndAtomIdx(),
#                                                   bnd.GetBondType())
        nr_bonds += bnd.BondType
        
    return nr_bonds


def main():
    aa_list=[a.lower() for a in ['Phe','Leu','Ile','Met','Val','Ser','Pro','Thr','Ala','Tyr', 'His', 'Gln', 'Asn', 'Lys', 'Asp', 'Glu', 'Cys', 'Trp', 'Arg', 'Gly']]
    
    matrix=''
    for a1_i,a1 in enumerate(aa_list,1):
        matrixrow=[]    
        for a2 in aa_list[:a1_i-1]:
            print os.getcwd()
            try:
                with open('{0}_{1}_mcs.out'.format(a1,a2)) as f: pass
            except IOError:
                subprocess.call('smsd -Q MOL -q amino/{0}.mol -T MOL -t amino/{1}.mol -o {0}_{1}.mol -g'.format(a1,a2),shell=True)
                subprocess.call('cp mcs.out {0}_{1}_mcs.out'.format(a1,a2),shell=True)
            
            f = open('{0}_{1}_mcs.out'.format(a1,a2), 'r')
            # first two lines should be like:
            #AtomContainer 1=    amino/phe.mol
            #AtomContainer 2=    amino/gly.mol
            line1=f.readline()
            line2=f.readline()
            if line1[-8:-5]!=a1 or line2[-8:-5]!=a2:
                raise
            size_sub=int(f.readline()[-3:])
            match1=[]
            match2=[]
            for i in arange(size_sub):
                # find the two numbers in the line and
                [mat1,mat2]=re.findall('\d+',f.readline())
                # append the first to match1, the second to match2, converted into integer
                match1.append(int(mat1))
                match2.append(int(mat2))
                
            
            # get mol structures without sanitation and removing H atoms
            m1=Molecule('amino/{0}.mol'.format(a1))
            m2=Molecule('amino/{0}.mol'.format(a2))            
            
            nr_bonds =  bonds_notin(m1,match1)
            nr_bonds += bonds_notin(m2,match2)
            
            print 'total nr of bonds between {0} and {1}: {2}'.format(a1,a2,nr_bonds)
            matrixrow.append('{:2.0f}'.format(nr_bonds))
    
        matrixrow.append(' 0')  # distance from itself is 0
        # latex formatting:
        matrix += '& '.join(matrixrow)  
        matrix += '\\\\ \n'
        
    print matrix
        

if __name__ == "__main__":
    main()
