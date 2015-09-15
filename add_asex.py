import re
import argparse
import sys
from numpy.random import poisson, randint

parser = argparse.ArgumentParser()

parser.add_argument('filename', type = str, help = 'filename containing ms-format loci (- for STDIN)')
parser.add_argument('numgens', type = int, help = 'number of generations of asex')
parser.add_argument('mutrate', type = float, help = 'per-base per-generation mutation probability')

args = parser.parse_args()

fin = open(args.filename, 'r') if args.filename != '-' else sys.stdin

possibleTypes = ['01', '10']

for line in fin:
    if re.match('@begin', line):
        segNumMuts = int(line.split()[1]) - 1

        line = fin.next()
        numBases = int(line.strip())
        
        mutPositions, mutTypes = [], []

        for i in xrange(segNumMuts):
            line = fin.next()
            spline = line.split()
            mutPosition, mutType = int(spline[0]), spline[1]

            mutPositions.append(mutPosition)
            mutTypes.append(mutType)

        poisMutRate = numBases * args.numgens * args.mutrate
        numAsexMuts = poisson(poisMutRate)

        asexPositions = [randint(0, numBases) for i in xrange(numAsexMuts)]
        asexTypes = [possibleTypes[randint(0,2)] for i in xrange(numAsexMuts)]

        mutPositions.extend(asexPositions)
        mutTypes.extend(asexTypes)

        mutPositions, mutTypes = zip(*sorted(zip(mutPositions, mutTypes)))

        newNumMuts = len(mutPositions)

        print "@begin {}".format(newNumMuts)
        print numBases
        for pos, typ in zip(mutPositions, mutTypes):
            print "{}\t{}".format(pos, typ)
    else:
        print line.strip()
