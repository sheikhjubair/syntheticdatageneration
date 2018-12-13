# -*- coding: utf-8 -*-
"""
Created on Sun Dec  9 13:51:25 2018

@author: Jubair
"""
import numpy as np
import ancestorDescendantGenerator as ad


#    
#    
duplicationProb = 1
inversionProb = 0
translocationProb = 0
transpositionProb = 0
adg = ad.ancestordescendantGenerator(2000,4, 1, 150, 5)
adg.getGenomes(duplicationProb, inversionProb, translocationProb, transpositionProb,3)
ancestralGenomes = adg.ancestralGenome
descendants = adg.descendants


filename = 'data/mduplication_600.txt'
genomeNumber = 1
with open(filename,'a') as file:
    file.write('# number of TD: {}\n'.format(adg.tdCount))
    file.write('# number of FD: {}\n'.format(adg.fdCount))         
    file.write('# number of Duplication: {}\n'.format(adg.duplicationCount))  
    file.write('# number of Inversion: {}\n'.format(adg.inversionCount))
    file.write('# number of Translocation: {}\n'.format(adg.translocationCount))
    file.write('# number of Transposition: {}\n'.format(adg.transpositionCount))
    file.write('\n')
    file.write('\n')
    file.write('\n')
#    file.write('Genome ' + str(genomeNumber) + ':\n')
#    
#    chromNumber = 1
#    for ancestorChrom in ancestralGenomes:
#        file.write('A_Chr_' + str(chromNumber)+ ' ')
#        np.savetxt(file, ancestorChrom, delimiter = ' ', fmt='%s', newline= ' ')
#        
#        file.write('\n')
#        chromNumber+=1
#    file.write('\n')
#    file.write('\n')
#    file.write('\n')
#    file.write('\n')
    genomeNumber+=1
    
    for ds in descendants:
        file.write('Genome ' + str(genomeNumber) + ':\n')
        genomeNumber+=1
        file.write('\n')
        chromNumber = 1
        for descendant in ds:
            file.write('D_Chr_' + str(chromNumber)+ ' ')
            np.savetxt(file, descendant, delimiter = ' ', fmt='%s', newline= ' ')
            file.write('\n')
            chromNumber+=1
        file.write('\n')
        file.write('\n')
        file.write('\n')
        file.write('\n')
        
        
        
        
        

#geneFamily = generateGeneFamily(1000)
#ancestChroms, linAncestChrom = generateAncestralGenome(geneFamily,5)
#descendants = generateDuplicatedChromosomesForGenome(linAncestChrom,1,250, 0.5, 1)

    