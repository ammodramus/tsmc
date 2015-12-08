import fileinput
from collections import Counter

for line in fileinput.input():
    line = line.strip()
    if line == '//':
        print line
        continue
    count = Counter()
    spline = [float(num) for num in line.split(',')]
    t2 = max(spline)
    for num in spline:
        if count[num] == 1:
            t3 = num
            break
        count[num] += 1
    t3 *= 2  # ms scaled by 4N
    t2 *= 2
    print "{},{}".format(t3,t2)
