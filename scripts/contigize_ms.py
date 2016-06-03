#!/usr/bin/env python
import re
import argparse
import sys
import numpy as np
import numpy.random as npr

def print_psmc_fasta_seq(positions, seqLen, blockLen, width, seqIdx):
    blocksAreHet = [False for i in xrange(seqLen/blockLen+1)]
    for pos in positions:
        block = pos / blockLen
        blocksAreHet[block] = True
    print ">{}_{}".format(seqIdx, seqLen)
    for i, blockIsHet in enumerate(blocksAreHet):
        if i > 0 and i % width == 0:
            print
        sys.stdout.write('K') if blockIsHet else sys.stdout.write('T'),
    print


parser = argparse.ArgumentParser()

parser.add_argument('ms', type = str, help = 'filename containing ms-format loci (- for STDIN)')
parser.add_argument('contiglengths', type = str, help = 'list of contig lengths to mimic')
parser.add_argument('--total-len', '-t', type = int, help = 'total length of genome to sample (bp)')

args = parser.parse_args()

contigLengths = []
with open(args.contiglengths, 'r') as conLenIn:
    for line in conLenIn:
        cLen = int(line.strip())
        contigLengths.append(cLen)
contigLengths = np.array(contigLengths, dtype = np.int)

fin = open(args.ms, 'r') if args.ms != '-' else sys.stdin

# ms-formatted file must have original ms command line with -t and -r flags
firstline = fin.next().strip()
spfirstline = firstline.split(' ')

numChroms = int(spfirstline[2])
chromLen = int(spfirstline[spfirstline.index('-r')+2])
theta = float(spfirstline[spfirstline.index('-t')+1])
genomeLength = numChroms * chromLen

if args.total_len is None:
    coveredLength = contigLengths.sum()
else:
    coveredLength = args.total_len

noncoveredLength = genomeLength-coveredLength

# get simulated contig lengths
sampledLengths = []
currentlyCovered = 0
while currentlyCovered < coveredLength:
    contigLen = npr.choice(contigLengths, size = 1, replace = True)[0]
    if currentlyCovered + contigLen > coveredLength:
        contigLen = coveredLength - currentlyCovered
    currentlyCovered += contigLen
    sampledLengths.append(contigLen)
sampledLengths = np.array(sampledLengths, dtype = np.int)

# break genome into uniformly spaced contigs

spaceStarts = np.hstack((np.array([0]),np.sort(npr.randint(0, noncoveredLength, size = sampledLengths.size))))
spaceLengths = np.diff(spaceStarts)

startPositions = np.cumsum(spaceLengths+sampledLengths)-sampledLengths
endPositions = startPositions+sampledLengths-1  # inclusive

# now get sequences!

# output needs to be psmcFA
curChrom = 0
allPositions = np.empty(0, dtype = np.int)
for line in fin:
    if re.match('^positions:', line):
        chromPositions = np.round(np.array([float(el) for el in line.split(' ')[1:]])*chromLen).astype(np.int) + chromLen*curChrom
        curChrom += 1
        allPositions = np.concatenate((allPositions, chromPositions))

# for each contig:
#   loop through positions until find one that isn't in range
#   then go to next contig

startEnds = zip(startPositions, endPositions)
startEndsNew = []
for i in xrange(len(startEnds)):
    start = startEnds[i][0]
    end = startEnds[i][1]
    if (start / chromLen) != (end / chromLen):
        start0 = start
        end0 = chromLen*((start / chromLen)+1)-1
        start1 = chromLen*((start / chromLen)+1)
        end1 = end
        startEndsNew.extend(((start0,end0), (start1,end1)))
    else:
        startEndsNew.append((start,end))

i = 0
curPolyIdx = 0
curPolyPos = allPositions[0]

while i < len(startEndsNew):
    start = startEndsNew[i][0]
    end = startEndsNew[i][1]
    contigPolys = []
    while curPolyPos <= end:
        if curPolyPos >= start and curPolyPos <= end:
            contigPolys.append(curPolyPos)
        curPolyIdx += 1
        try:
            curPolyPos = allPositions[curPolyIdx]
        except IndexError:
            break # no more positions, print
    contigLen = end - start + 1
    contigRelPolys = [pos-start for pos in contigPolys]
    print_psmc_fasta_seq(contigRelPolys, contigLen, 100, 60, i)
    i += 1

