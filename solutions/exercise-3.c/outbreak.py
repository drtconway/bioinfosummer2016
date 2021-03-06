from pykmer.basics import kmers
from pykmer.file import readFasta

import sys
import gzip

def jaccard(xs, ys):
    xz = len(xs)
    yz = len(ys)
    i = 0
    j = 0
    d = 0
    b = 0
    while i < xz and j < yz:
        x = xs[i]
        y = ys[j]
        if x < y:
            d += 1
            i += 1
            continue
        if x > y:
            d += 1
            j += 1
            continue
        b += 1
        i += 1
        j += 1
    d += xz - i
    d += yz - j
    return float(b) / float(b + d)

K = 27

kmerSets = []
for fn in sys.argv[1:]:
    print >> sys.stderr, fn
    with gzip.open(fn) as f:
        xs = set([])
        for (nm, seq) in readFasta(f):
            xs |= set(kmers(K, seq, True))
        xs = list(xs)
        xs.sort()
        kmerSets.append((fn, xs))

for i in range(len(kmerSets)):
    for j in range(i + 1, len(kmerSets)):
        d = jaccard(kmerSets[i][1], kmerSets[j][1])
        print '%s\t%s\t%g' % (kmerSets[i][0], kmerSets[j][0], d)
