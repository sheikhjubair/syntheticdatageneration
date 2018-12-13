# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 00:32:08 2018

@author: Jubair
"""

import numpy as np

class ancestordescendantGenerator:
    def __init__(self, geneFamilyLength, chromInAncestor, dupliLength, numEvolutionEvents, seed):
        self.geneFamilyLength = geneFamilyLength
        self.geneFamily = np.chararray((self.geneFamilyLength,), itemsize =10, unicode=True)
        self.seed = seed
        self.dupliLength = dupliLength
        self.numEvolutionEvents = numEvolutionEvents
        self.inversionCount = 0
        self.duplicationCount = 0
        self.ancestralGenome = []
        self.descendants = []
        self.chromInAncestor = chromInAncestor
        self.translocationCount = 0
        self.transpositionCount = 0
        self.tdCount = 0
        self.fdCount = 0
        
        
        
    def generateGeneFamily(self):
        ind = 0
        for _i in range(0, self.geneFamilyLength):
            gene = 'g'+ str(_i+1)
            self.geneFamily[ind] = gene
            ind+=1
    
    
    
    def generateAncestralGenome(self, numChromosomes):
        minGenesPerChrom = np.int(self.geneFamily.shape[0] / numChromosomes)
        
        ancestralChroms = []
        currentGeneIndex = 0
        linearAncestorChrom = self.createReverseOrientation(self.geneFamily)
        chromType = np.chararray((numChromosomes,), itemsize =10, unicode=True)
        for chromNum in range(numChromosomes):
            self.seed = self.seed + 10
            np.random.seed(self.seed)
            if numChromosomes != chromNum + 1:
                chrom = linearAncestorChrom[currentGeneIndex:currentGeneIndex + minGenesPerChrom]
            else:
                chrom = linearAncestorChrom[currentGeneIndex:]
                
            currentGeneIndex = currentGeneIndex + minGenesPerChrom
            
            isCircular = np.random.choice(2)
        
            if isCircular:
                chromType[chromNum] = ')'
            else:
                chromType[chromNum] = '|'
            ancestralChroms.append(chrom)
            
        return ancestralChroms, chromType
        

    def createTranslocation(self, chrom1, chrom2, maxTranslocateLength):
        #maxTranslocateLength must be greater than 2
        self.seed = self.seed + 10
        np.random.seed(self.seed)
        chrom1StartPoint = np.random.randint(0, chrom1.shape[0]-maxTranslocateLength -1)
        chrom1EndPoint = chrom1StartPoint + np.random.randint(2, maxTranslocateLength)
        chrom1TranslocateGenes = chrom1[chrom1StartPoint: chrom1EndPoint]
        
        chrom2StartPoint = np.random.randint(0, chrom2.shape[0]-maxTranslocateLength -1)
        chrom2EndPoint = chrom2StartPoint + np.random.randint(2, maxTranslocateLength)
        chrom2TranslocateGenes = chrom2[chrom2StartPoint: chrom2EndPoint]
        
        chrom1= np.delete(chrom1, np.arange(chrom1StartPoint, chrom1EndPoint))
        chrom2 = np.delete(chrom2, np.arange(chrom2StartPoint, chrom2EndPoint))
        chrom1 = np.insert(chrom1, chrom1StartPoint, chrom2TranslocateGenes)
        chrom2 = np.insert(chrom2, chrom2StartPoint, chrom1TranslocateGenes)
        self.translocationCount +=1
            
        return chrom1, chrom2
    
    def createTransposition(self, chromosome, maxTranspositionLength):
        self.seed = self.seed + 10
        np.random.seed(self.seed)
        startPos = np.random.randint(0, chromosome.shape[0] - maxTranspositionLength)
        endPos = startPos + np.random.randint(1, maxTranspositionLength)
        transpositionGenes = chromosome[startPos: endPos]
        chromosome = np.delete(chromosome, np.arange(startPos, endPos))
        transpositionSite = np.random.randint(0, chromosome.shape[0])
        chromosome = np.insert(chromosome, transpositionSite, transpositionGenes)
        self.transpositionCount +=1
        return chromosome
        

    def createInversion(self, chromosome):
        self.seed = self.seed + 10
        np.random.seed(self.seed)
        startPos = np.random.randint(1, chromosome.shape[0] -2)
        endPos = np.random.randint(startPos+1, chromosome.shape[0])
        
        genesForInversions = chromosome[startPos: endPos]
        genesForInversions = np.flip(genesForInversions, 0)
        ind = 0
        invertedGenes = np.chararray((genesForInversions.shape[0],), itemsize =10, unicode=True)
        for gene in genesForInversions:
            if gene[0] == '-':
                invertedGenes[ind] = gene[1:]
            else:
                invertedGenes[ind] = '-'+gene
            ind+=1
        chromosome[startPos: endPos] = invertedGenes
        
        
        self.inversionCount+=1
        
        return chromosome
    
    
    def createDuplication(self, chromosome, dupliLength):
        self.seed = self.seed+ 10
        np.random.seed(self.seed)
        duplicationPos = np.random.randint(0, chromosome.shape[0] - dupliLength +1)
        isTandemDuplication = np.random.randint(0,2)
        duplicatedGenes = chromosome[duplicationPos: duplicationPos + dupliLength]
        
        if isTandemDuplication:
            chromosome = np.insert(chromosome, duplicationPos + dupliLength, duplicatedGenes, axis = 0)
            self.tdCount +=1
        else:
            breakPointPosition = np.random.randint(duplicationPos+ dupliLength, chromosome.shape[0]+1)
            chromosome = np.insert(chromosome, breakPointPosition, duplicatedGenes, axis = 0)
            self.fdCount += 1
        self.duplicationCount+=1
        return chromosome

    def createReverseOrientation(self,chrom):
        numRevOrientation = np.int(chrom.shape[0] *0.25)
        revInd = np.random.choice(np.arange(0, chrom.shape[0]), numRevOrientation)
        chrom[revInd] = '-'+ chrom[revInd]
        
        return chrom
    
    def getGenomes(self, dupliProb, inverProb, transloProb, tanspoProb, numDescendants):
        self.generateGeneFamily()
        self.ancestralGenome, chromType = self.generateAncestralGenome(self.chromInAncestor)
        
        
        for i in range(numDescendants):
            descendant = list.copy(self.ancestralGenome)
            self.seed = self.seed + 25 + i
            for evolEvent in range(self.numEvolutionEvents):
                self.seed = self.seed + 10
                np.random.seed(self.seed)
                eventType = np.random.choice(4, p= [dupliProb, inverProb, transloProb, tanspoProb])
                chromosomeNumbers = np.random.choice(self.chromInAncestor, 2, replace = False)
                if eventType == 0:
                    descendant[chromosomeNumbers[0]] = self.createDuplication(descendant[chromosomeNumbers[0]], self.dupliLength)
                elif eventType == 1:
                    descendant[chromosomeNumbers[0]] = self.createInversion(descendant[chromosomeNumbers[0]])
                elif eventType == 2:
                    descendant[chromosomeNumbers[0]], descendant[chromosomeNumbers[1]] = self.createTranslocation(descendant[chromosomeNumbers[0]], descendant[chromosomeNumbers[1]], 10)
                else:
                    descendant[chromosomeNumbers[0]] = self.createTransposition(descendant[chromosomeNumbers[0]],10)
            
            for chromNum in range(len(descendant)):
                descendant[chromNum] = np.append(descendant[chromNum], chromType[chromNum])
            self.descendants.append(descendant)
            
        for chromNum in range(len(self.ancestralGenome)):
            self.ancestralGenome[chromNum] = np.append(self.ancestralGenome[chromNum], chromType[chromNum])
                    
        

    
        
        
    