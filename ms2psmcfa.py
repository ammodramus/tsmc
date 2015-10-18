import fileinput
import sys

BLOCK_LEN = 100
OUTPUT_WIDTH = 60

positions = []
chromCount = 1
for line in fileinput.input():
    if fileinput.isfirstline():
        splitcmd = line.strip().split(' ')
        numSegSites = int(splitcmd[splitcmd.index('-r')+2])
    line = line.strip()
    if line[:11] == 'positions: ':
        chromPositions = [float(el) for el in line.split(' ')[1:]]
        chromPositions = [int(numSegSites*pos) for pos in chromPositions]
        blocksAreHet = [False for i in xrange(numSegSites/BLOCK_LEN+1)]
        for chromPos in chromPositions:
            block = chromPos / BLOCK_LEN
            blocksAreHet[block] = True
        print '>{}'.format(chromCount)
        for i, blockIsHet in enumerate(blocksAreHet):
            if i > 0 and i % OUTPUT_WIDTH == 0:
                print
            sys.stdout.write('K') if blockIsHet else sys.stdout.write('T'),
        print
