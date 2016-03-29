#!/usr/bin/env python
import sys
import subprocess as sp
import numpy as np
import numpy.random as npr

usage = '''
Usage: diptripgenome.py NUMCHROMS TRIPTIME DIPTIME LAMBDAD MSOPTIONS...
Simulates an asexual triploid genome derived from a fertilized diploid asexual
genome.

Options:
    numchroms  number of chromosomes to simulate
    triptime   amount of time between present and fertilization of diploid asexual
    diptime    amount of time between fertilization of diploid asexual and
               origin of diploid asexual
    lambdad    population size during the diploid asexual part of the history,
               relative to N0
    msoptions  additional options to ms, including -r, -t, and any demographic history
               (-r and -t must be provided, and the demographic history is for
               the sexual ancestor population only)

All times are measured in units of 4*N0, as in ms. scrm must be available in
the system's PATH. Does not support multiple populations.
'''

numRequiredOptions = 4

if '-h' in sys.argv or '--help' in sys.argv:
    print usage
    sys.exit(0)

if len(sys.argv) < numRequiredOptions+1:
    print '\nToo few arguments'
    print usage
    sys.exit(1)

if '-r' not in sys.argv:
    print '\nMust provide -r argument, formatted as in ms\n'
    print usage
    sys.exit(1)

if '-t' not in sys.argv:
    print '\nMust provide -t argument\n'
    print usage
    sys.exit(1)

numChroms = int(sys.argv[1])
triptime = float(sys.argv[2])
diptime = float(sys.argv[3])
lamd = float(sys.argv[4])

# process msoptions
msoptions = sys.argv[(numRequiredOptions+1):]
popnSizeChangeIdxs = [i for i, x in enumerate(msoptions) if x in ['-eN', '-eG']]

for popIdx in popnSizeChangeIdxs:
    oldtime = float(msoptions[popIdx+1])
    newtime = oldtime + diptime
    msoptions[popIdx+1] = str(newtime)

try:
    rho = float(msoptions[msoptions.index('-r')+1])
    chromLen = int(msoptions[msoptions.index('-r')+2])
    theta = float(msoptions[msoptions.index('-t')+1])
except:
    print '\nBad ms options'
    print usage
    sys.exit(1)

scrmCmd = [str(el) for el in ['scrm', '3', numChroms, '-eI', diptime, 2,
    '-eI', '0.0', '1']]
scrmCmd.extend([str(el) for el in msoptions])

out = str(sp.check_output(scrmCmd)).split('\n')

outIter = iter(out)

for line in outIter:
    line = line.strip()
    if line == '//':
        line = outIter.next().strip() # segsites
        numSegsites = int(line.split(' ')[1])
        line = outIter.next().strip() # positions: 0.0419644 0.043589 ...
        positions = np.array([float(el) for el in line.split(' ')[1:]])
        line = outIter.next().strip() # 01011110010100010110101001000001001...
        firstHap = np.array([int(el) for el in list(line)])
        line = outIter.next().strip()
        secondHap = np.array([int(el) for el in list(line)])
        line = outIter.next().strip()
        thirdHap = np.array([int(el) for el in list(line)])

        
        # asex mutations for two lineages that were asexual first
        numAsexMuts01 = npr.poisson(2.0*(diptime+triptime) * 2.0*theta/2.0)  # note 2.0 multiplier
        numAsexMuts2 = npr.poisson(2.0*triptime * theta/2.0)
        numAsexMuts = numAsexMuts01 + numAsexMuts2

        asexMutPositions01 = npr.rand(numAsexMuts01)  # npr.rand(n) simulates n Unif(0,1)'s
        asexMutPositions2 = npr.rand(numAsexMuts2)

        asexTypes01 = npr.choice(np.array([0,1]), size=numAsexMuts01, replace = True) # which lineage gets the asex mutation

        asexFirstHapTypes = np.zeros(numAsexMuts01, dtype = np.int)
        asexSecondHapTypes = np.zeros(numAsexMuts01, dtype = np.int)
        asexThirdHapTypes = np.zeros(numAsexMuts01, dtype = np.int)
        asexFirstHapTypes[asexTypes01 == 0] = 1
        asexSecondHapTypes[asexTypes01 == 1] = 1

        asexFirstHapTypes = np.concatenate((asexFirstHapTypes, np.zeros(numAsexMuts2, dtype = np.int)))
        asexSecondHapTypes = np.concatenate((asexSecondHapTypes, np.zeros(numAsexMuts2, dtype = np.int)))
        asexThirdHapTypes = np.concatenate((asexThirdHapTypes, np.ones(numAsexMuts2, dtype = np.int)))

        positions = np.concatenate((positions, asexMutPositions01, asexMutPositions2))
        firstHap = np.concatenate((firstHap, asexFirstHapTypes))
        secondHap = np.concatenate((secondHap, asexSecondHapTypes))
        thirdHap = np.concatenate((thirdHap, asexThirdHapTypes))

        sortedIdxs = positions.argsort()
        positions = positions[sortedIdxs]
        firstHap = firstHap[sortedIdxs]
        secondHap = secondHap[sortedIdxs]
        thirdHap = thirdHap[sortedIdxs]
        newNumSegSites = positions.size

        print '//'
        print 'segsites: {}'.format(newNumSegSites)
        print 'positions: {}'.format(' '.join([str(el) for el in positions]))
        print ''.join([str(el) for el in list(firstHap)])
        print ''.join([str(el) for el in list(secondHap)])
        print ''.join([str(el) for el in list(thirdHap)])

    else:
        print line
