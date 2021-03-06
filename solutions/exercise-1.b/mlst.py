from pykmer.basics import  kmers
from pykmer.file import readFasta

import gzip
import sys

K = 27

# Step 1
# Index the alleles
idx = {}
lens = {}
for loc in sys.argv[1:-1]:
    with open(loc) as f:
        for (nm, seq) in readFasta(f):
            xs = set(kmers(K, seq, True))
            for x in xs:
                if x not in idx:
                    idx[x] = set([])
                idx[x].add(nm)
            lens[nm] = len(xs)
                
# Because k-mers may occur multiple times in a sequence,
# we keep track of those in the index that we've already
# seen, so we don't reprocess them.
#
seen = set([])

# This dictionary keeps all the partially computed results:
# alleles for which we've seen at least one k-mer, which is
# to say, we record the number of k-mers we need to find to
# complete the allele.
#
partial = {}

results = []
with gzip.open(sys.argv[-1]) as f:
    for (_, seq) in readFasta(f):
        for x in kmers(K, seq, True):
            if x not in idx:
                continue

            # k-mer was in the index,
            # but we've already processed it.
            if x in seen:
                continue
            seen.add(x)

            for nm in idx[x]:
                if nm not in partial:
                    # First k-mer for the allele.
                    # Initialize with the number
                    # of k-mers in the allele.
                    partial[nm] = lens[nm]

                partial[nm] -= 1

                if partial[nm] == 0:
                    # Bingo! We found all the
                    # k-mers in the allele.
                    results.append(nm)

results.sort()
print '\t'.join(results)
