from pykmer.basics import kmers
from pykmer.file import readFasta

import sys
import gzip

def jaccard(xs, ys):
    return float(len(xs & ys)) / float(len(xs | ys))

K = 27

kmerSets = []
for fn in sys.argv[1:]:
    print >> sys.stderr, fn
    with gzip.open(fn) as f:
        xs = set([])
        for (nm, seq) in readFasta(f):
            xs |= set(kmers(K, seq, True))
        print len(xs)
        kmerSets.append((fn, xs))

for i in range(len(kmerSets)):
    for j in range(i + 1, len(kmerSets)):
        d = jaccard(kmerSets[i][1], kmerSets[j][1])
        print '%s\t%s\t%g' % (kmerSets[i][0], kmerSets[j][0], d)
