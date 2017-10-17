import sys
sys.path.append("../lib/")

from re import split

class AlignError(Exception):
    def __init__(self, msg = ""):
        self.msg = str(msg)
    def __str__(self):
        return self.msg

def readScoringMatrix(scoreMatrixFile):
    """Read in a GenBank-formated scoring matrix and
    return a scoring matrix in dictionary form"""
    with open(scoreMatrixFile) as fp:
        line = fp.readline()
        while line[0] == '#':
            line = fp.readline()

        char_list = split("\s+", line.strip())
        D = {}

        for line in fp:
            arr = split("\s+", line.strip())
            D[arr[0]] = {char_list[i]:float(val)  for i,val in enumerate(arr[1:])}
    return D

def scoreAlignment(aln1, aln2, S, gap):
    """Compute the score of an alignment relative to scoring matrix S and gap
    penelty gap (represented as a non-negative number) -- for testing."""

    if {type(aln1),type(aln2)} != {str}:
        raise AlignError("Aligned sequences are not strings")
    if len(aln1) != len(aln2):
        raise AlignError("Aligned sequences are not of the same length")
    if not (set(aln1) | set(aln2) <= set(S.keys()) | {'-'}):
        raise AlignError("Aligned sequences contain bad characters")
    
    return sum([-gap if '-' in {x,y} else S[x][y]  for x,y in zip(aln1,aln2)])
