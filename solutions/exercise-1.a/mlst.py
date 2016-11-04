from pykmer.basics import  kmers
from pykmer.file import readFasta

import gzip
import sys

K = 27
alleles = []
for loc in sys.argv[1:-1]:
    with open(loc) as f:
        for (nm, seq) in readFasta(f):
            alleles.append((nm, set(kmers(K, seq, True))))

xs = set([])
with gzip.open(sys.argv[-1]) as f:
    for (nm, seq) in readFasta(f):
        xs |= set(kmers(K, seq, True))

results = []
for (nm, ys) in alleles:
    if len(xs & ys) == len(ys):
        results.append(nm)

results.sort()
print '\t'.join(results)
