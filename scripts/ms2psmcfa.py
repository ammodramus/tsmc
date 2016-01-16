#!/usr/bin/env python
import fileinput
import numpy as np
import sys

BLOCK_LEN = 100
OUTPUT_WIDTH = 60

positions = []
chromCount = 1
inp = fileinput.input()
for line in inp:
    if fileinput.isfirstline():
        splitcmd = line.strip().split(' ')
        numSites = int(splitcmd[splitcmd.index('-r')+2])
        numHaps = int(splitcmd[1])
        continue
    line = line.strip()
    if line[:11] == 'positions: ':
        if numHaps == 2:   # diploid
            chromPositions = [float(el) for el in line.split(' ')[1:]]
            chromPositions = [int(numSites*pos) for pos in chromPositions]
            blocksAreHet = [False for i in xrange(numSites/BLOCK_LEN+1)]
            for chromPos in chromPositions:
                block = chromPos / BLOCK_LEN
                blocksAreHet[block] = True
            print '>{}'.format(chromCount)
            for i, blockIsHet in enumerate(blocksAreHet):
                if i > 0 and i % OUTPUT_WIDTH == 0:
                    print
                sys.stdout.write('K') if blockIsHet else sys.stdout.write('T'),
            print
            chromCount += 1
        elif numHaps == 3:   # triploid
            chromPositions = np.array([float(el) for el in line.split(' ')[1:]], dtype = np.float64)
            chromPositions = (numSites * chromPositions).astype(np.int)
            chromSNPtypes = np.zeros(len(chromPositions), dtype = np.int)
            for i in xrange(3):
                line = inp.readline().strip()
                chromSNPtypes += np.array([int(el) for el in line], dtype = np.int)
            numBlocks = numSites/BLOCK_LEN+1
            blockBits = np.zeros((numBlocks, 2))
            for SNPpos, SNPtype in zip(list(chromPositions), list(chromSNPtypes)):
                SNPblock = SNPpos / BLOCK_LEN
                blockBits[SNPblock, SNPtype-1] = 1
            blockTypes = np.zeros(numBlocks)
            blockTypes += blockBits[:,0]
            blockTypes += 2*blockBits[:,1]
            print '>{}'.format(chromCount)
            for i, blockType in enumerate(blockTypes.astype(np.int)):
                if i > 0 and i % OUTPUT_WIDTH == 0:
                    print
                sys.stdout.write(str(blockType)),
            print
            chromCount += 1


                
