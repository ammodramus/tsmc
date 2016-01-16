#!/usr/bin/env python
import re
import argparse
import sys
import numpy as np
import numpy.random as npr

parser = argparse.ArgumentParser()

parser.add_argument('filename', type = str, help = 'filename containing ms-format loci (- for STDIN)')
parser.add_argument('divtime', type = float, help = 'number of generations since transition from sexual reproduction, in units of 2*N0')
parser.add_argument('--triploid', action='store_true')

args = parser.parse_args()

fin = open(args.filename, 'r') if args.filename != '-' else sys.stdin

# theta is 4*N0*mu in ms
# here theta is multiplied by time in units of 2*N0
#

# ms-formatted file must have original ms command line with -t and -r flags
firstline = fin.next().strip()
spfirstline = firstline.split(' ')
numSites = int(spfirstline[spfirstline.index('-r')+2])
theta = float(spfirstline[spfirstline.index('-t')+1])
if not args.triploid and spfirstline[1] != '2':
    raise ValueError('expected 2 haplotypes in ms output format')
if args.triploid and spfirstline[1] != '3':
    raise ValueError('expected 3 haplotypes in ms output format for triploids')

print firstline
if not args.triploid:
    for line in fin:
        line = line.strip()
        if re.match('//', line):
            line = fin.next().strip()  # segsites: 12398
            segNumMuts = int(line.split(' ')[1])
            line = fin.next().strip()  # positions: 0.0005 0.0011 
            positions = np.array([float(el) for el in line.split(' ')[1:]], dtype = np.float64)
            line = fin.next().strip() # 01011110010100010110101001000001001...
            firstHap = np.array([int(el) for el in list(line)])
            line = fin.next().strip()
            # secondHap = line

            # simulate asexual / Meselson-effect mutations
            numAsexMuts = npr.poisson(args.divtime * theta)
            asexMutPositions = npr.rand(numAsexMuts)  # npr.rand(n) simulates n Unif(0,1)'s
            asexTypes = npr.choice(np.array([0,1]), size=numAsexMuts, replace = True)

            # adding the asexual mutations
            positions = np.concatenate((positions, asexMutPositions))
            firstHap = np.concatenate((firstHap, asexTypes))
            sortedIdxs = positions.argsort()
            positions = positions[sortedIdxs]
            firstHap = firstHap[sortedIdxs]
            secondHap = np.bitwise_xor(firstHap, 1)
            newNumSegSites = positions.size

            # print the new locus
            print '//'
            print 'segsites: {}'.format(newNumSegSites)
            print 'positions: {}'.format(' '.join([str(el) for el in positions]))
            print ''.join([str(el) for el in list(firstHap)])
            print ''.join([str(el) for el in list(secondHap)])

        else:
            print line

elif args.triploid:
    for line in fin:
        line = line.strip()
        if re.match('//', line):
            line = fin.next().strip()  # segsites: 12398
            segNumMuts = int(line.split(' ')[1])
            line = fin.next().strip()  # positions: 0.0005 0.0011 
            positions = np.array([float(el) for el in line.split(' ')[1:]], dtype = np.float64)
            line = fin.next().strip() # 01011110010100010110101001000001001...
            firstHap = np.array([int(el) for el in list(line)])
            line = fin.next().strip()
            secondHap = np.array([int(el) for el in list(line)])
            line = fin.next().strip()
            thirdHap = np.array([int(el) for el in list(line)])

            # simulate asexual / Meselson-effect mutations
            numAsexMuts = npr.poisson(args.divtime * 3.0*theta/2.0)  # note 3/2 multiplier
            asexMutPositions = npr.rand(numAsexMuts)  # npr.rand(n) simulates n Unif(0,1)'s
            asexTypes = npr.choice(np.array([0,1,2]), size=numAsexMuts, replace = True)
            asexFirstHapTypes = np.zeros(asexTypes.shape[0], dtype = np.int)
            asexSecondHapTypes = np.zeros(asexTypes.shape[0], dtype = np.int)
            asexThirdHapTypes = np.zeros(asexTypes.shape[0], dtype = np.int)
            asexFirstHapTypes[asexTypes == 0] = 1
            asexSecondHapTypes[asexTypes == 1] = 1
            asexThirdHapTypes[asexTypes == 2] = 1


            # adding the asexual mutations
            positions = np.concatenate((positions, asexMutPositions))
            firstHap = np.concatenate((firstHap, asexFirstHapTypes))
            secondHap = np.concatenate((secondHap, asexSecondHapTypes))
            thirdHap = np.concatenate((thirdHap, asexThirdHapTypes))

            sortedIdxs = positions.argsort()

            positions = positions[sortedIdxs]
            firstHap = firstHap[sortedIdxs]
            secondHap = secondHap[sortedIdxs]
            thirdHap = thirdHap[sortedIdxs]
            newNumSegSites = positions.size

            # print the new locus
            print '//'
            print 'segsites: {}'.format(newNumSegSites)
            print 'positions: {}'.format(' '.join([str(el) for el in positions]))
            print ''.join([str(el) for el in list(firstHap)])
            print ''.join([str(el) for el in list(secondHap)])
            print ''.join([str(el) for el in list(thirdHap)])

        else:
            print line
