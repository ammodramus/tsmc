#!/usr/bin/env python
import re
import argparse
import sys
import numpy as np
import numpy.random as npr

parser = argparse.ArgumentParser()

parser.add_argument('filename', type = str, help = 'filename containing ms-format loci (- for STDIN)')
parser.add_argument('divtime', type = float, help = 'number of generations since transition from sexual reproduction, in units of 2*N0')

args = parser.parse_args()

fin = open(args.filename, 'r') if args.filename != '-' else sys.stdin

firstline = fin.next().strip()
spfirstline = firstline.split(' ')
numSites = int(spfirstline[spfirstline.index('-r')+2])
theta = float(spfirstline[spfirstline.index('-t')+1])
if spfirstline[1] != '2':
    raise ValueError('expected 2 haplotypes in ms output format')

possibleTypes = ['01', '10']

print firstline
for line in fin:
    line = line.strip()
    if re.match('//', line):
        line = fin.next().strip()
        segNumMuts = int(line.split(' ')[1])
        line = fin.next().strip()
        positions = np.array([float(el) for el in line.split(' ')[1:]], dtype = np.float64)
        line = fin.next().strip()
        firstHap = [int(el) for el in list(line)]
        line = fin.next().strip()
        secondHap = line
        numAsexMuts = npr.poisson(args.divtime * theta)
        asexMutPositions = npr.rand(numAsexMuts)
        asexTypes = npr.choice(np.array([0,1]), size=numAsexMuts, replace = True)
        positions = np.concatenate((positions, asexMutPositions))
        firstHap = np.concatenate((firstHap, asexTypes))
        sortedIdxs = positions.argsort()
        positions = positions[sortedIdxs]
        firstHap = firstHap[sortedIdxs]
        secondHap = np.bitwise_xor(firstHap, 1)
        newNumSegSites = positions.size
        print '//'
        print 'segsites: {}'.format(newNumSegSites)
        print 'positions: {}'.format(' '.join([str(el) for el in positions]))
        print ''.join([str(el) for el in list(firstHap)])
        print ''.join([str(el) for el in list(secondHap)])
    else:
        print line

