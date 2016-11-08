from pykmer.basics import  kmers
from pykmer.file import readFasta
import pykmer.kfset as kfset

import gzip
import sys

def isFasta(nm):
    """does this filename look like a FASTA file?"""
    if nm.endswith(".fa"):
        return True
    if nm.endswith(".fas"):
        return True
    if nm.endswith(".fasta"):
        return True
    if nm.endswith(".fna"):
        return True
    return False

K = 27

# Step 1
# Index the alleles
idx = {}
lens = {}
profArg = 1
for loc in sys.argv[1:]:
    if not isFasta(loc):
        break
    with open(loc) as f:
        for (nm, seq) in readFasta(f):
            xs = set(kmers(K, seq, True))
            for x in xs:
                if x not in idx:
                    idx[x] = set([])
                idx[x].add(nm)
            lens[nm] = len(xs)
    profArg += 1

# Step 2,
# Index the profile tuples.
profiles = {}
headers = None
with open(sys.argv[profArg]) as f:
    for l in f:
        t = l.strip().split('\t')
        if headers is None:
            headers = t
            continue
        st = t[0]
        aroC = "AROC" + t[2]
        dnaN = "DNAN" + t[3]
        hemD = "HEMD" + t[4]
        hisD = "HISD" + t[5]
        purE = "PURE" + t[6]
        sucA = "SUCA" + t[7]
        thrA = "THRA" + t[8]
        prof = (aroC, dnaN, hemD, hisD, purE, sucA, thrA)
        profiles[prof] = st

for fn in sys.argv[profArg+1:]:
    seen = set([])
    partial = {}
    results = []
    (_, xs) = kfset.read(fn)
    for (x,_) in xs:
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
    r = tuple(results)
    if r in profiles:
        print fn + '\t' + profiles[r]
    else:
        print fn + '\t' + 'unknown:' + '\t'.join(results)
