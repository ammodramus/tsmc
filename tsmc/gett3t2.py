import fileinput
from collections import Counter

for line in fileinput.input():
    count = Counter()
    line = line.strip()
    spline = [float(num) for num in line.split(',')]
    t2 = max(spline)
    for num in spline:
        if count[num] == 1:
            t3 = num
            break
        count[num] += 1
    print "{},{}".format(t3,t2)
