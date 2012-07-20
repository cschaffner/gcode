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
        if self.EndAtomIdx == None or self.BondType<0 or self.BondType>2 or self.BondType == None:
            raise

class Molecule:
    def __init__(self,molfile):
        """ imports bonds from a molfile """
        self.name=molfile[-7:-4]
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
    
    def GetBondByAtomIdx(self,BeginAtomIdx,EndAtomIdx):
        for b in self.Bonds:
            if b.BeginAtomIdx==BeginAtomIdx and b.EndAtomIdx==EndAtomIdx:
                return b
            elif b.BeginAtomIdx==EndAtomIdx and b.EndAtomIdx==BeginAtomIdx:
                return b
        return False
    
    def GetNrBonds(self):
        nr_bonds=0
        for bnd in self.Bonds:
            nr_bonds += bnd.BondType
        return nr_bonds 

def bonds_notin(m1,match1,m2,match2):    
    """ 
    takes as input molecules m1,m2 and list of atom-indices match1 and match2
    where the maximal common subgraph of m1 and m2 is isomorphic to
    the atoms in match1 of m1
    and the atoms in match2 of m2
        
    returns the number of bonds that are in molecule m1, 
    but are not part of the maximal common subgraph
    """
    
    nr_bonds=0
    for bnd in m1.GetBonds():
        # check if that bond is possibly part of the maximal common subgraph
        if bnd.BeginAtomIdx in match1 and bnd.EndAtomIdx in match1:
            # check if that bond also exists in the other molecule m2
            # (it has to if it should be part of the maximal common subgraph)
            m2Begin=match2[match1.index(bnd.BeginAtomIdx)]
            m2End=match2[match1.index(bnd.EndAtomIdx)]
            m2Bond=m2.GetBondByAtomIdx(m2Begin,m2End)
            if m2Bond == False:
                 # if the bond does not exist in the other molecule m2,
                 # it cannot be part of the maximal common subgraph, so we have to add it to the bond_count
                 print 'bond does not exist in m2, so not part of mcs: adding extra {0}'.format(bnd.BondType)
                 nr_bonds += bnd.BondType
            elif bnd.BondType == 2 and m2Bond.BondType == 1:
                # otherwise (i.e. if the bond exists in m2)
                # but here it is a double-bond, and in m2 it's only a single-bond, we have to
                # add one to the bond count
                print 'simple bond in subgraph, but double in real: adding an extra bond' 
                nr_bonds += 1
            # if the bond is part of the maximal common subgraph, we do not have to count it
            continue
        else: 
            # if the bond is NOT part of the maximal common subgraph, we have to count it anyway
            nr_bonds += bnd.BondType

    print 'bond difference from {0} to mcs of {0} and {1}: {2}'.format(m1.name,m2.name,nr_bonds)
    return nr_bonds


def main():
    aa_list=[a.lower() for a in ['Phe','Leu','Ile','Met','Val','Ser','Pro','Thr','Ala','Tyr', 'His', 'Gln', 'Asn', 'Lys', 'Asp', 'Glu', 'Cys', 'Trp', 'Arg', 'Gly']]
    
    matrix=''
    for a1_i,a1 in enumerate(aa_list,1):
        matrixrow=[]    
        for a2 in aa_list[:a1_i-1]:
            try:
                with open('{0}_{1}_mcs.out'.format(a1,a2)) as f: pass
            except IOError:
                subprocess.call('smsd -Q MOL -q amino/{0}.mol -T MOL -t amino/{1}.mol -o {0}_{1}.mol -g -m'.format(a1,a2),shell=True)
                subprocess.call('cp mcs.out {0}_{1}_mcs.out'.format(a1,a2),shell=True)
            
            # read in structures
            m1=Molecule('amino/{0}.mol'.format(a1))
            m2=Molecule('amino/{0}.mol'.format(a2))
                                        
            f = open('{0}_{1}_mcs.out'.format(a1,a2), 'r')
            # first two lines should be like:
            #AtomContainer 1=    amino/phe.mol
            #AtomContainer 2=    amino/gly.mol
            line1=f.readline()
            line2=f.readline()
            if line1[-8:-5]!=a1 or line2[-8:-5]!=a2:
                raise
            size_sub=int(f.readline()[-3:])

            # go through all possible solutions of matchings and find the one
            # with the minimal number of bonds to change to go from one to the other            
            min_bonds=m1.GetNrBonds()+m2.GetNrBonds()
            while True:
                # first an empty line
                firstline=f.readline()
                if firstline == '':
                    # end of file reached
                    break
                line=f.readline()
                if line[:8] != 'Solution':
                    print 'incorrect format of mcs file'
                    raise                
                else:
                    print 'checking solution: {0}'.format(line[-3:])
                match1=[]
                match2=[]
                for i in arange(size_sub):
                    # find the two numbers in the line and
                    [mat1,mat2]=re.findall('\d+',f.readline())
                    # append the first to match1, the second to match2, converted into integer
                    match1.append(int(mat1))
                    match2.append(int(mat2))                                                                        
                nr_bonds =  bonds_notin(m1,match1,m2,match2)
                nr_bonds += bonds_notin(m2,match2,m1,match1)
                if nr_bonds < min_bonds:
                    min_bonds = nr_bonds
                # another empty line 
                f.readline()                   
                # then the "//" line
                f.readline()
            
            print 'total nr of bonds between {0} and {1}: {2}\n'.format(a1,a2,nr_bonds)
            matrixrow.append('{:2.0f}'.format(nr_bonds))
    
        matrixrow.append(' 0')  # distance from itself is 0
        # latex formatting:
        matrix += '& '.join(matrixrow)  
        matrix += '\\\\ \n'
        
    print matrix
        

if __name__ == "__main__":
    main()
