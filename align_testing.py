from align_util import *
import align
import re

def test_alignment(a1, a2, s1, s2, score, correct_score, S, g):
    assert re.sub("-", "", a1) == s1, "First aligned string is not based on the original string"
    assert re.sub("-", "", a2) == s2, "Second aligned string is not based on the original string"
    assert score == correct_score, "The returned score is not optimal"
    assert score == scoreAlignment(a1, a2, S, g), "The returned score is not the score of the returned alignment"


def test_NeedlemanWunsch():
    
    # Test 1
    print("Test 1:")
    s1 = "CCC"
    s2 = "CCC"
    S = readScoringMatrix("DNA1.txt")
    score, a1, a2 = align.NeedlemanWunsch(s1, s2, S, 1)
    test_alignment(a1, a2, s1, s2, score, 3, S, 1)
    print("Passed\n")

    # Test 2
    print("Test 2:")
    s1 = "CCC"
    s2 = "CTC"
    S = readScoringMatrix("DNA1.txt")
    score, a1, a2 = align.NeedlemanWunsch(s1, s2, S, 1)
    test_alignment(a1, a2, s1, s2, score, 2, S, 1)
    print("Passed\n")

    # Test 3
    print("Test 3:")
    s1 = "CTC"
    s2 = "CCC"
    S = readScoringMatrix("DNA1.txt")
    score, a1, a2 = align.NeedlemanWunsch(s1, s2, S, 1)
    test_alignment(a1, a2, s1, s2, score, 2, S, 1)
    print("Passed\n")

    # Test 4
    print("Test 4:")
    s1 = "CCC"
    s2 = "CC"
    S = readScoringMatrix("DNA1.txt")
    score, a1, a2 = align.NeedlemanWunsch(s1, s2, S, 1)
    test_alignment(a1, a2, s1, s2, score, 1, S, 1)
    print("Passed\n")
    
    # Test 5
    print("Test 5:")
    s1 = "CC"
    s2 = "CCC"
    S = readScoringMatrix("DNA1.txt")
    score, a1, a2 = align.NeedlemanWunsch(s1, s2, S, 1)
    test_alignment(a1, a2, s1, s2, score, 1, S, 1)
    print("Passed\n")

    # Test 6
    print("Test 6:")
    s1 = "AGTA"
    s2 = "ACA"
    S = readScoringMatrix("DNA2.txt")
    score, a1, a2 = align.NeedlemanWunsch(s1, s2, S, 2)
    test_alignment(a1, a2, s1, s2, score, 7, S, 2)
    print("Passed\n")

    # Test 7
    print("Test 7:")
    s1 = "ACA"
    s2 = "AGTA"
    S = readScoringMatrix("DNA2.txt")
    score, a1, a2 = align.NeedlemanWunsch(s1, s2, S, 2)
    test_alignment(a1, a2, s1, s2, score, 7, S, 2)
    print("Passed\n")

    # Test 8
    print("Test 8:")
    s1 = "WWWGGWWW"
    s2 = "WWWLWWW"
    S = readScoringMatrix("Blosum62.txt")
    score, a1, a2 = align.NeedlemanWunsch(s1, s2, S, 2)
    test_alignment(a1, a2, s1, s2, score, 60, S, 2)
    print("Passed\n")



    
    
    