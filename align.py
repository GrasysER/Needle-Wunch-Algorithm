import sys
sys.path.append("../lib/")

from align_testing import *
from align_util import *

def test():
    s1 = "ACA"
    s2 = "AGTA"
    S = readScoringMatrix("DNA2.txt")    
    return NeedlemanWunsch(s1, s2, S, 2)
    

def NeedlemanWunsch(seq1, seq2, S, gap):
    """
    Input:
    * seq1, seq2: Two *strings* to be aligned.  (These could be either dna sequences of protein sequences -- it shouldn't matter.)
    * S: a scoring dictionary such that for any characters a,b: S[a][b] is the miss-match penalty for aligning an a to a b.
      -- You may assume that every character contained in seq1 or seq2 is defined in S.
    * gap: The gap score, listed as a non-negative number.  (e.mg. If gap=4, we will deduct 4 points for each -.)
    Output: a tupple containing:
    * the optimal score (as a float)
    * the first aligned string
    * the second aligned string
    """
    gap = -1 * gap
    s1 = list(seq1)

    s2 = list(seq2)
 
    #Score
    n = len(s1)
    m = len(s2)
    M = [[0 for i in range(m+1)] for j in range(n+1)]
    
    for j in range(m+1):
        M[0][j] = -j
    for i in range(1, n+1):
        M[i][0] = -i 
        for j in range(1, m+1):
            t1 = M[i-1][j-1] + S[s1[i-1]][s2[j-1]]
            t2 = M[i-1][j] + gap
            t3 = M[i][j-1] + gap
            M[i][j] = max(t1, t2, t3)
        
    
    
    #Backtrace
    A = B = "" 
    i = len(s1) 
    j = len(s2) 
    print("Length of strings:", i, j)

    while (i > 0 and j > 0):
            Score = M[i][j]
            diag = M[i - 1][j - 1]
            up = M[i][j - 1]
            left = M[i - 1][j]
            if (Score == diag + S[s1[i-1]][s2[j-1]]):
                    A = s1[i-1] + A
                    B = s2[j-1] + B
                    i -= 1
                    j -= 1
            elif (Score == left + gap):
                    A = s1[i-1] + A
                    B = "-" + B
                    i -= 1
            elif (Score == up + gap):
                    A = "-" + A
                    B = s2[j-1] + B
                    j -= 1
         
    while (i > 0):
            A = s1[i] + A
            B = "-" + B
            i -= 1
    while (j > 0):
            A = "-" + A
            B = s2[j] + B
            j -= 1
    
    
    return (M[len(s1)][len(s2)], A, B)    
    
   