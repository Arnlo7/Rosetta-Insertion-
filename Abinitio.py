#Aaron Lao
#02-251 Great Ideas in Computational Biology
import os
from tkinter import *
from tkinter.filedialog import askopenfilename
import copy
import math
import random
import string
import numpy
import vpython


#Helper Functions
def collapseWhitespace(s):
    result=""
    inWhiteSpace=False
    for c in s:
        if(c==" "or c=="\n" or c=="\t"):
            if(not inWhiteSpace):
                result+=" "
                inWhiteSpace=True
        else:
            result+=c
            inWhiteSpace=False
    return result

def distance(coords1,coords2):
    term1=(coords2[0]-coords1[0])**2
    term2=(coords2[1]-coords1[1])**2
    term3=(coords2[2]-coords2[2])**2
    return math.sqrt(term1+term2+term3)

##Object Classes###
class Residue(object):
    def __init__(self,residue,psi,phi,omega):
        #In order, N-C-C-O
        self.atomList=[]
        self.residue=residue
        self.psi=psi
        self.phi=phi
        self.omega=omega
        self.replaced=False
        self.centroid= AminoAtom(residue,f"A{residue}",0,0,0)

    def addAtom(self,atom):
        self.atomList.append(atom)

    #should only be called after alpha carbon updated
    def updateCentroid(self):
        self.centroid.x=self.atomList[1].x
        self.centroid.y=self.atomList[1].y+3
        self.centroid.z=self.atomList[1].z

class AminoAtom(object):
    def __init__(self,residue,element, x, y, z):
        self.vradiusDict=dict()
        self.vradiusDict["H"]=.12
        self.vradiusDict["O"]=.152
        self.vradiusDict["C"]=.17
        self.vradiusDict["Ca"]=.17
        self.vradiusDict["N"]=.155
        #Vander Waals Radius: approximate, sum vander waals radii of heavy
        #components. VERY approximate
        #Alanine: CH3
        self.vradiusDict["AA"]= 3*.12+.17
        #Arginine: C4 H10 N2
        self.vradiusDict["AR"]=4*.17+10*.12+2*.155
        #Asparagine
        self.vradiusDict["AN"]=2*.17+.155+.17+.152
        #Aspartic Acid
        self.vradiusDict["AD"]=2*.17+2*.12+2*.152
        #Cysteine
        self.vradiusDict["AC"]=3*.12+.17
        #Glutamic Acid
        self.vradiusDict["AE"]=3*.17+4*.12+2*.152
        #Glutamine
        self.vradiusDict["AQ"]=3*.17+6*.12+.155+.152
        #Glycine
        self.vradiusDict["AG"]=.12
        #Histidine
        self.vradiusDict["AH"]=5*.17+2*.155
        #Isoleucine
        self.vradiusDict["AI"]=4*.17+9*.12
        #Leucine
        self.vradiusDict["AL"]=4*.17+9*.12
        #Lysine
        self.vradiusDict["AK"]=4*.17+11*.12+.155
        #Methionine
        self.vradiusDict["AM"]=3*.17+7*.12
        #Phenylalanine
        self.vradiusDict["AF"]=7*.17+7*.12
        #Proline
        self.vradiusDict["AP"]=3*.17+6*.12
        #Serine
        self.vradiusDict["AS"]=.17+3*.12+.152
        #Thereonine
        self.vradiusDict["AT"]=2*.17+5*.12+.152
        #Tryptophan
        self.vradiusDict["AW"]=10*.17+7*.12+.155
        #Tyrosine
        self.vradiusDict["AY"]=7*.17+7*.12+.152
        #Valine
        self.vradiusDict["AV"]=3*.17+7*.12
        self.residue=residue
        self.element=element
        self.x=x
        self.y=y
        self.z=z
        self.vradius=self.vradiusDict[self.element]

class CountObject(object):
    def __init__(self):
        self.count=0

#class ThreeNeighbor(object):
    #def __init__(self,):

class NineNeighbor(object):
    def __init__(self,insertion):
        self.insertion=insertion
        self.neighbors=[]

    def addNeighbor(self,neighbor):
        self.neighbors.append(neighbor)

class Fragment(object):
    def __init__(self,residue,phi,psi,omega):
        "Letter representing AA"
        self.residue=residue
        self.phi=phi
        self.psi=psi
        self.omega=omega

class ExtractFragments(object):
    def __init__(self,toggle):
        if(toggle==0):

            pass
        elif(toggle==1):
            self.data="""
"""
        else:
            self.txtInput()

    def build9Data(self):
        splitList= self.data.split("position")
        currResidue=1
        neighborList=[]
        for stringData in splitList:
            stringData=stringData.strip()
            nestList= stringData.splitlines()
            tempString=""
            for i in range(len(nestList)):
                if(repr(nestList[i])!=repr("") and repr(nestList[i]!=repr(" "))):
                    if(i!=0):
                        tempString+=f"{nestList[i]}\n"
            tempList=(tempString.strip()).splitlines()
            positionList=[]
            for i in range(0,len(tempList),9):
                nine = NineNeighbor(currResidue)
                for j in range(i,i+9):
                    currString=tempList[j].strip()
                    residue=currString[14]
                    angleString=currString[18:44]
                    angleList=(collapseWhitespace(angleString).strip()).split(" ")
                    phi=float(angleList[0])
                    psi=float(angleList[1])
                    omega=float(angleList[2])
                    frag = Fragment(residue,phi,psi,omega)
                    nine.addNeighbor(frag)
                positionList.append(nine)
            if(positionList!=[]):
                currResidue+=1
                neighborList.append(positionList)
        return neighborList
    #If the fasta is n long, neighbor list has n entries
    #Each entry in neighbor list is a list with 200 9Neighbors

    def txtInput(self):
        filename=askopenfilename()
        self.file=open(filename)
        self.data=(self.file.read()).strip()
        self.file.close()

class SequenceReader(object):
    #Assume that the fragments were run on the same protein
    def __init__(self,toggle,neighborList):
        self.neighborList=neighborList
        #In terms of angstroms
        #Numbers obtained from ref 1
        self.lengthDict=dict()
        self.lengthDict["N-C"]=1.329
        #self.lengthDict["C-N1"]=1.341 #if proline
        self.lengthDict["C-O"]=1.231
        self.lengthDict["Ca-C"]=1.525
        #self.lengthDict["Ca-C1"]=1.516 #if gly
        #Not sure what beta carbons are, leave out for now
        #self.lengthDict["Ca-Cb2"]=1.53
        #self.lengthDict["Ca-Cb1"]=1.521 #If Ala
        #self.lengthDict["Ca-Cb2"]=1.540 #If Ile, Thr, Val
        self.lengthDict["N-Ca"]=1.458 
        #self.lengthDict["N-Ca1"]=1.451 #if gly
        #self.lengthDict["N-Ca2"]=1.466 #if pro
        #define O-N as the same as C-N bond
        self.lengthDict["O-N"]=1.329

        #Bond Angles
        self.angleDict=dict()
        #Assuming regardless of side chain configuration, 
        #and with Ca as some carbon bonded to a hydrogen
        #in tetrahedral arrangement (from ref #1)
        self.angleDict["N-Ca-C"]= 111.2
        #self.angleDict["Ca-N-C"]= 123.0
        #self.angleDict["N-C-Ca"]= 114.0
        #self.angleDict["N-C-O"]= 125.0
        self.angleDict["Ca-C-O"]=121.0
        self.angleDict["C-N-Ca"]=123.0
        #just define C-O-N as equivalent to  O-C-N, just a checking purpose
        #self.angleDict["O-C-N"]=125.0
        self.angleDict["C-O-N"]=125.0

        
        self.atomOrder=["N","Ca","C","O"]
        self.currIndex=0
        self.atomList=[]
        if(toggle==0):

            pass
        elif(toggle==1):
            self.data="""
"""
        else:
            self.txtInput()
        #Assume all torsion angles are 180 (for extended conformation)


    def run(self):
        #Initialization Step:
        #self.extendedConformation()
        #self.monteCarlo()
        #selected=self.select9frag()
        #for neighbor in selected[0].neighbors:
            #print(neighbor.psi)
        #self.insertNeighbor(self.residues,selected[0],selected[1])
        return self.outputPymol(self.fragmentSampler()[1])

###SCORING FUNCTION####
    def coarseRosettaScore(self):
        energy=0
        #ermList=[]
        #Residue Environment (solvation)
        #env
        #termList.append(env)
        #Residue Pair Interactions(electrostatics, disulfides)
        #pair
        #termList.append(pair)
        #Strand Pairing (Hydrogen Bonding)
        #SS
        #termList.append(SS)
        #Strand Arrangement into Sheets
        #sheet
        #termList.append(sheet)
        #Helix Strand Packing
        #HS
        #termList.append(HS)
        #Radius of Gyration
        rg=self.gyrationRadius()
        #termList.append(rg)
        #CB density
        #cbeta
        #termList.append(cbeta)
        #Steric Repulsion
        vdw=self.stericRepulsion()
        #termList.append(vdw)        
        #for i in range(len(weightVector)):
            #energy+=weightVector[i]*termList[i]
        return rg+vdw
        #return energy

# The vdw term represents only steric repulsion and not attractive van der Waals forces
# (those are modeled in terms rewarding compact structures, such as the rg term; local 
# interactions are implicitly included from fragments). It is calculated over pairs of 
# atoms only in cases where: 1. the interatomic distance is less than the sum of the atoms' 
# van der Waals radii, and 2. the interatomic distance does not depend on the torsion angles 
# of a single residue. (Rohl et al. 2004)
    def stericRepulsion(self):
        total=0
        for res1 in (self.residues):
            for res2 in (self.residues):   
                r=res1.centroid.vradius+res2.centroid.vradius
                d=distance((res1.centroid.x,res1.centroid.y,res1.centroid.z),(res2.centroid.x,res2.centroid.y,res2.centroid.z))
                total+=(((r**2)-(d**2))**2)/r
        return total

#Radius ofGyration:
    def gyrationRadius(self):
        total=0
        for res1 in (self.residues):
            for res2 in (self.residues):
                d=distance((res1.centroid.x,res1.centroid.y,res1.centroid.z),(res2.centroid.x,res2.centroid.y,res2.centroid.z))
                total+=math.sqrt(d**2)
        return total

##MONTE CARLO ANNEALING SEARCH
    #General Idea:
    #1. Start with protein in fully extended conformation (not sure how i'll represent it for now)
    #2. Randomly select a 9 residue fragment insertion window.
    #3. Select random top 25 fragments in list. Replace torsion angles
    #4. Retain moves that decrease energy, use metropolis criterion on others
    def fragmentSampler(self):
        bestScore= math.inf
        #atomList is now in extended conformation
        self.extendedConformation()
        replacedCount= CountObject()
        bestConfiguration=self.residues
        #while(replacedCount.count<len(self.residues)-100)
        for i in range(100):
            nineFragment, insertion=self.select9frag()
            selected=self.select9frag()
            oldConformation=copy.deepcopy(self.residues)
            self.insertNeighbor(self.residues,selected[0],selected[1],replacedCount)
            #Conformation Updated
            #Stage 1: Continues until all initial torsion angles have been replaced
            score=self.coarseRosettaScore()
            #Only consider steric overlap term
            #elif(_<2000):
                #score=rosettaFunction(self.residues)
            #Consider all terms except cbeta, ra
            #elif(_<22000):
                #score=rosettaFunction(self.residues)
            #Consider strand pairing
            #elif(_<28000):
                #score=rosettaFunction(self.residues)
            if(score<bestScore*2):
                bestScore=score
            else:
                self.residues=oldConformation
            #else:
                #if(metropolisCriterion(score, other parameters)):
                    #bestScore=score
                    #currentConformation=childConformation
        return bestScore, self.residues
            
    def select9frag(self):
        #returns selected fragment and its insertion spot
        #arbitrarily position 2 and sequence 1
        randomLocation=random.randint(1,len(self.residues)-13)
        random25=random.randint(0,24)
        return self.neighborList[randomLocation][random25],randomLocation
    #def select9frag(self):
        #ninemerList=[]
        #for i in range(self.queryLength-9+1):
            #self.neighborList.append(self.query[i:i+9])
        #rand=random.randint(0,len(ninemerList)-1)
        #return ninemerList[rand], rand

###INITIALIZATION: brings protein to native structure
    def extendInitialize(self,firstResidue,resObject):
        #In the extended conformation, all torsion angles are 180
        #Starts with N-terminus
        first="N"
        nitrogenAtom=AminoAtom(firstResidue,first,0,0,0)
        self.atomList.append(nitrogenAtom)
        self.currIndex+=1
        #Arbitrarily set initial angle as zero
        second="Ca"
        carbonAtom=AminoAtom(firstResidue,second,self.lengthDict[f"{first}-{second}"],0,0)
        self.atomList.append(carbonAtom)
        self.currIndex+=1
        resObject.addAtom(nitrogenAtom)
        resObject.addAtom(carbonAtom)
        #returns the third atom as the latest atom
        return nitrogenAtom, carbonAtom

    def extendCycle(self,residue,prev2, prev1,resObject):
        prev1El=prev1.element
        prev2El=prev2.element
        oldX=prev1.x
        oldY=prev1.y
        if(prev2El=="O" and prev1El=="N"):
            new="Ca"
            angle=self.angleDict["C-N-Ca"]
            length=self.lengthDict["N-C"]
            newX=oldX+length*math.cos(180-angle)
            newY=oldY+length*math.sin(180-angle)
            newAtom=AminoAtom(residue,new,newX,newY,0)
            self.atomList.append(newAtom)
            resObject.addAtom(newAtom)
            self.currIndex+=1
        else:
            new=self.atomOrder[self.currIndex%4]
            angle=self.angleDict[f"{prev2El}-{prev1El}-{new}"]
            newX=oldX+self.lengthDict[f"{prev1El}-{new}"]*math.cos(math.radians(180-angle))
            newY=oldY+self.lengthDict[f"{prev1El}-{new}"]*math.sin(math.radians(180-angle))
            newAtom=AminoAtom(residue,new,newX,newY,0)
            self.atomList.append(newAtom)
            resObject.addAtom(newAtom)
            self.currIndex+=1
        return prev1, newAtom

    def extendResidue(self,residue,prev2,prev1,newRes):
        for _ in range(4):
            prev2,prev1=self.extendCycle(residue,prev2,prev1,newRes)
        return prev2,prev1

    def extendedConformation(self):
        #for each residue, initialize nitrogen, acarbon, carbon, oxygen 
        self.residues=[]
        newRes= Residue(self.data[0],-180,-180,-180)
        prev2,prev1=self.extendInitialize(self.data[0],newRes)
        prev2,prev1=self.extendCycle(self.data[0],prev2,prev1,newRes)
        prev2,prev1=self.extendCycle(self.data[0],prev2,prev1,newRes)
        self.residues.append(newRes)
        for residue in self.data[1:]:
            if(residue!="\n"):
                newRes= Residue(residue,-180,-180,-180)
                prev2,prev1=self.extendResidue(residue,prev2,prev1,newRes)
                self.residues.append(newRes)
        return self.atomList
  
###TORSION ANGLE INSERTION####
    ##Only apply these functions after computing 2D conformation!
    def translateResidues(self,Residues,translationTuple):
        for residue in Residues:
            for atom in residue.atomList:
                atom.x=atom.x-translationTuple[0]
                atom.y=atom.y-translationTuple[1]
                atom.z=atom.z-translationTuple[2]
    #This is not a void function, returns new coordinates
    #This is after 0,0,0 already in place
    def singleRotation(self,coordinates1,coordinates2,angle):
        v1= vpython.vector(coordinates2[0],coordinates2[1],coordinates2[2])
        axis = vpython.vector(coordinates1[0],coordinates1[1],coordinates1[2])
        result= vpython.rotate(v1,math.radians(angle),axis)
        return (result.x,result.y,result.z)

    #This is a void function, changes values in Residues
    #currAtom = [N->0,Ca->1,C->2]
    #currResidue 
    #If Omega rotation, currAtom will be N (curr1Coor)
    #If Phi rotation, currAtom will be Ca
    #IF Psi rotation, currAtom will be C
    def rotateConformation(self,Residues,currResIndex,currAtom,prevCoor,curr1Coor,angle):
        self.translateResidues(Residues,curr1Coor)
        #Applied on the remaining atoms in the residue
        #Will be applied on Ca, C, O
        coordinates1=prevCoor
        for i in range(3-currAtom):
            atom=Residues[currResIndex].atomList[currAtom+1+i]
            coordinates2=(atom.x,atom.y,atom.z)
            atom.x,atom.y,atom.z=self.singleRotation(coordinates1,coordinates2,180-angle)
        Residues[currResIndex].updateCentroid()
        #Applied on all remaining residues in the chain
        for i in range(currResIndex+1,len(Residues)):
            for atom in Residues[i].atomList:
                coordinates2=(atom.x,atom.y,atom.z)
                atom.x,atom.y,atom.z=self.singleRotation(coordinates1,coordinates2,180-angle)
            Residues[i].updateCentroid()

    #Starts with nitrogen
    def insertFragment(self,Residues,currResidue,previousRes,fragment,location,countOb):
        #Update angles:
        currResidue.omega=fragment.omega
        currResidue.phi=fragment.phi
        currResidue.psi=fragment.psi
        if(currResidue.replaced):
            currResidue.replaced=True
            countOb.count+=1
        prevC=previousRes.atomList[2]
        prevCCoor=(prevC.x,prevC.y,prevC.z)
        currN=currResidue.atomList[0]
        currNCoor=(currN.x,currN.y,currN.z)
        currCa=currResidue.atomList[1]
        currCaCoor=(currCa.x,currCa.y,currCa.z)
        currC=currResidue.atomList[2]
        currCCoor=(currC.x,currC.y,currC.z)
        nextN=Residues[location+1].atomList[0]
        nextNCoor=(nextN.x,nextN.y,nextN.z)
        #Omega Parameters: 
        #currAtom=0, prevCoor=previous carbon,curr1Coor= nitrogen, angle= previous omega
        #Buggy must fix
        self.rotateConformation(Residues,location,0,prevCCoor,currCaCoor,fragment.omega)
        prevC=previousRes.atomList[2]
        prevCCoor=(prevC.x,prevC.y,prevC.z)
        currN=currResidue.atomList[0]
        currNCoor=(currN.x,currN.y,currN.z)
        currCa=currResidue.atomList[1]
        currCaCoor=(currCa.x,currCa.y,currCa.z)
        currC=currResidue.atomList[2]
        currCCoor=(currC.x,currC.y,currC.z)
        nextN=Residues[location+1].atomList[0]
        nextNCoor=(nextN.x,nextN.y,nextN.z)
        #Phi Parameters
        #currAtom=1, prevCoor=current Nitrogen, curr1Coor= current Ca, angle= curr phi
        self.rotateConformation(Residues,location,1,currNCoor,currCCoor,fragment.phi)
        prevC=previousRes.atomList[2]
        prevCCoor=(prevC.x,prevC.y,prevC.z)
        currN=currResidue.atomList[0]
        currNCoor=(currN.x,currN.y,currN.z)
        currCa=currResidue.atomList[1]
        currCaCoor=(currCa.x,currCa.y,currCa.z)
        currC=currResidue.atomList[2]
        currCCoor=(currC.x,currC.y,currC.z)
        nextN=Residues[location+1].atomList[0]
        nextNCoor=(nextN.x,nextN.y,nextN.z)
        #Psi Parameters
        #currAtom=2, prevCoor=current Ca, curr1Coor= current C, angle= curr psi
        self.rotateConformation(Residues,location,2,currCaCoor,nextNCoor,fragment.psi)

    #INPUT: A list Residues of residue objects, a fragment, and an insertion point
    #Insertion point specifies which amino acid
    def insertNeighbor(self,Residues,nineNeighbor,insertion,countOb):
        insertion=insertion-1
        if(insertion==0 or insertion==len(Residues)-1): #edge case
            pass
        else:
            previousRes= Residues[insertion-1]
            for i in range(insertion,insertion+9):
                currRes=Residues[i]
                self.insertFragment(Residues,currRes,previousRes,nineNeighbor.neighbors[i-insertion],i,countOb)
                previousRes=currRes

###OUTPUT ANSWER###
    def txtInput(self):
        filename=askopenfilename()
        self.file=open(filename)
        self.data=self.file.read().strip()
        self.data.strip()
        self.file.close()

    def outputPymol(self,residues):
        newString=""
        count=1
        for residue in residues:
            for atom in residue.atomList:
                resNum=count%4+1
                displayX="%0.3f" % (atom.x)
                displayY="%0.3f" % (atom.y)
                displayZ="%0.3f" % (atom.z)
                if(atom.element=="Ca"):
                    displayElem="C"
                else:
                    displayElem=atom.element
                newString+=f"ATOM      {count}  {displayElem}   {atom.residue} A  {resNum}      {displayX}  {displayY}  {displayZ}  1.00 50  {displayElem} 0"
                count+=1
                newString+="\n"
        return newString
###TESTER###
print("Testing...")
sys.setrecursionlimit=3000
#Select fragment file
fragment = ExtractFragments(2)
#select fasta sequence
reader = SequenceReader(2,fragment.build9Data())
answerString = reader.run()
answerFile= open("answer.pdb","w+")
answerFile.write("")
answerFile.write(f"""
{answerString.strip()}
""")
#answerFile.close()
print("Done!")




"""
##Amino Acid Data Dictionary
#class AminoData(object):
    #def __init__(self):

#Fragment Sampler:
#Library of possible fragments obtained from Rosetta Database
#Standard cycle size is 2000 cycles
#class FragmentSampler(object):

    #Residue3 and Residue9 are the Rosetta outputted fragment matches
    #Both have the form of a very long string
    #def __init__(self):
        #print("Select 3-residue fragments and 9-residue fragments")
        #self.resInput()
        #print("Select query sequence")
        #self.queryInput()
        #self.queryLength=len(self.query)
        #Convert to some sort of list

    #def resInput(self):
        #residue3
        #filename=askopenfilename()
        #res3File=open(filename)
        #self.res3=res3File.read()
        #res3File.close()
        #residue9
        #filename=askopenfilename()
        #res9File=open(filename)
        #self.res9=res9File.read()
        #res9File.close()

    #def queryInput(self):
        #filename=askopenfilename()
        #queryFile=open(filename)
        #QueryFile should be txt, obtain from FASTA format
        #self.query=queryFile.read()
        #queryFile.close()
    

    #def select9mer(self):
        #ninemerList=[]
        #n-k+1 possible kmers
        #for i in range(self.queryLength-9+1):
            #nineMerList.append(self.query[i:i+9])
        #rand=random.randint(0,len(ninemerList)-1)
        #return ninemerList[rand], rand
    #Takes an integer database type
    #def selectFragment(self,database):
        #if(database==3):

        #elif(database==9):
            #tempList=[]
            #for i in range(25):
                #Take highest 25 scores from the fragment database
                #tempList.append(self.res9[i])
            #rand=random.randint(0,24)
            #return tempList[rand]
    
    #Takes a fragment and startPosition. Recalculates the
    #torsion angles at this position with the new torsion angles
    #of the fragment
    #def replaceTorsion(self,fragment,startPos):
    





#Scoring Function:
#Torsion space representation:
#Protein backbone specified by torsion angles, modifications
#occur within the torsion space. 
#Cartesian space protein reprsentation generated with
#atomic coordinates for all heavy atoms in protein backbone
#(assuming idieal bond lengths and angles for individual residues)

#Reduced Energy Function:
#1. Each side chain represented by centroid located at side-chain 
#center of masss. Centroid determined by averaging over observed side-chain 
#conformations in known protein structures. Atomic coordinates for ALL side chain atoms
#are utilized. (will have to be external, uses protein database)
#2. Do not explicitly represent solvent
#3. Assume all bond lengths and bond angles are fixed
#4. Represent protein backbone using torsion angnles
    #Computes score of current conformation according to many parameters
    #Rough grained search, follows above assumptions
    #def rosettaFunction(many parameters lmao)


#References:
#1. Engh R A & Huber R (1991). Accurate bond and angle parameters for X-ray protein structure refinement. Acta Cryst., A47, 392-400
"""
