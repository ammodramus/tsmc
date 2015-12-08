import numpy as np
import fileinput
import sys

if len(sys.argv) != 3:
    raise Exception("Usage: python tallyt3t2.py timegraphy input")

def get_time_index(time):
    global timepoints
    for idx in range(len(timepoints)):
        if timepoints[idx+1] > time:
            return idx

def get_rowcol_index(i, j):
    global n
    idx = i*n-i*(i-1)/2+j;
    return idx;

timeFilename = sys.argv[1]
inp = sys.argv[2]

if inp == '-':
    inp = sys.stdin
else:
    inp = open(inp, 'r')

timepoints = []
timeInp = open(timeFilename, 'r')
for line in timeInp:
    line = line.strip()
    time = float(line.split('\t')[1])
    timepoints.append(time)

timepoints.append(float('inf'))

n = len(timepoints)-2 

t3s = []
t2s = []
for line in inp:
    line = line.strip()
    if line == '//':
        t3s.append([])
        t2s.append([])
        continue
    t3, t2 = [float(el) for el in line.split(',')]
    t3s[-1].append(t3)
    t2s[-1].append(t2)

numStates = (n+1)*(n+2)/2

counts = np.zeros(shape=(numStates,numStates), dtype = np.int)

for t3sChrom, t2sChrom in zip(t3s, t2s):
    s3 = t3sChrom[0]
    s2 = t2sChrom[0]
    i = get_time_index(s3)
    j = get_time_index(s2)

    for t3, t2 in zip(t3sChrom[1:], t2sChrom[1:]):
        if t3 == s3 and t2 == s2:
            continue
        k = get_time_index(t3)
        l = get_time_index(t2)
        if i == k and j == l:
            continue
        rowIdx = get_rowcol_index(i,j)
        colIdx = get_rowcol_index(k,l)
        counts[rowIdx,colIdx] += 1
        s3 = t3
        s2 = t2

probs = counts.astype(np.double)/counts.sum()
np.savetxt(sys.stdout, probs, delimiter = ',')
