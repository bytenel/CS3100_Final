import re
from dfa import *
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Alphabet
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable
from random import randint
from random import randrange

''' mk_dna_dfa():
    Takes a species number and obtains the start codons and a set of stop codons and generates a corresponding DFA
    Facts: 
    1) Non coding sequences are allowed  - A valid string contains a sequence of non coding and coding
        sequences (coding sequences are designated by start and stop codons)
    2) Must contain at least one coding sequence
    3) Does not allow for point or frame shift mutations - they can be added later
'''
def mk_dna_dfa(species):
    Q = {'S0', 'S0A', 'S0B', 'Smid', 'SmidA', 'SmidB', 'Send', 'SendA', 'SendB'}
    Sigma = {'A','T','G','C'}
    Delta = dict()
    F = {'Send', 'SendA', 'SendB'}
    q0 = 'S0'

    table = CodonTable.unambiguous_dna_by_id[species]
    start_codons = table.start_codons
    stop_codons = table.stop_codons

    count = 0

    """There are three loops that we must account for"""
    for n in range(3):

        if n==0:
            """First loop - Non Coding loop"""
            origin = 'S0'           #loop begins at 'S0' (start of DFA)
            destination = 'Smid'    #loop ends at 'Smid' (midpoint of DFA)
            sideA = 'S0A'           
            sideB = 'S0B'
            codons = start_codons   #loop will end once a start codon is read
        elif n==1:
            """Second loop - Coding Loop"""
            origin = 'Smid'         #loop begins at 'Smid' (midpoint of DFA)
            destination = 'Send'    #loop ends at 'Send' (endpoint of DFA)
            sideA = 'SmidA'         
            sideB = 'SmidB'
            codons = stop_codons    #loop will end once a stop codon is read
        elif n==2:
            """Third loop - Non Coding Loop"""
            origin = 'Send'        #loop begins at 'Send' (endpoint of DFA)
            destination = 'Smid'   #goes back 'Smid' (midpoint - start of coding sequence)
            sideA = 'SendA'
            sideB = 'SendB'
            codons = start_codons   #loop will end once a start codon is read

        paths = list()      #contains the codons that have been taken care of so far
        newstates1 = list() #first tier of new states (incorrect entries here will go to side state B)
        newstates2 = list() #second tier of new states (incorrect entries here will go to the origin state)
        
        for i in codons:
            """Make a path for all of the codons from the origin to the destination"""
            nucleotide1 = i[0]  #get the first nucleotide of the current codon
            nucleotide2 = i[1]  #get the second nucleotide of the current codon                       
            nucleotide3 = i[2]  #get the third nucleotide of the current codon

            """Check to see if we can merge paths with another codon"""
            solved = False
            for j in paths:                         #check all of the codons that have been taken care of already
                if not solved :       #dont check the current codon against itself!
                    result = get_difference(i, j)   #get the difference in nucleotides between the two codons
                    if result[0]==1:                #if they are off by one nucleotide, no new states - add a path somewhere
                        if result[1]==1 and not Delta.has_key((origin,nucleotide1)):
                            """If they differ by the first codon"""
                            state_to_join = Delta[(origin, j[0])]               #get the state lead to by the other codon's first nucleotide
                            Delta.update({(origin,nucleotide1): state_to_join}) #add this codon's first nucleotide to that definition
                            solved = True
                        elif result[2]==1:
                            """If they differ by the second codon"""
                            stateA = Delta[(origin, j[0])]                      #get the state lead to by the other codon's first nucleotide
                            state_to_join = Delta[(stateA, j[1])]               #get the second state in the path
                            Delta.update({(stateA,nucleotide2): state_to_join}) #add the different nucleotide to the path in between the first and second state
                            solved = True
                        elif result[3]==1:
                            """If they differ by the third codon"""
                            stateA = Delta[(origin, j[0])]                  #get the state lead to by the other codon's first nucleotide
                            stateB = Delta[(stateA, j[1])]                  #get the second state in the path
                            Delta.update({(stateB,nucleotide3): destination})    #add the different nucleotide to the path in between the second state and the destination
                            solved = True
            if not solved:  #if we are not able to merge this with another codon's path we must add at least one state
                for j in paths:
                    if i is not j and not solved:
                        result = get_difference(i, j)
                        if result[0]==2:    #if the codons are off by two (excluding cases where the 1st and last one are different)
                            if result[1]==1 and result[2]==1:
                                """If the first two nucleotides are different"""
                                stateA = Delta[(origin, j[0])]          #get the state lead to by the other codon's first nucleotide
                                state_to_join = Delta[(stateA, j[1])]   #get the last state in the path of the comparison codon

                                count += 1                              
                                new_state = 'S' + str(count)                        #new state name
                                Delta.update({(origin, nucleotide1): new_state})    #path from origin to new state with current nucleotide
                                Q.add(state1)                                       #must update Q
                                newstates1.append(state1)                           #this is in the first tier of states

                                Delta.update({(origin,nucleotide1): state_to_join})   #merge this path with the last state in the other path
                                solved = True

                            elif result[2]==1 and result[3]==1:
                                """if the last two nucleotides are different"""
                                state_to_join = Delta[(origin, j[0])]     #get the state that the first nucleotide leads to

                                count += 1                              
                                new_state = 'S' + str(count)                            #new state name
                                Delta.update({(state_to_join, nucleotide2): new_state}) #path from shared state to new state with current nucleotide
                                Q.add(state1)                                           #must update Q
                                newstates2.append(state1)                               #this is in the second tier of states

                                Delta.update({(new_state,nucleotide3): destination})   #path from this state to destination
                                solved = True
            if not solved:
                """If it is not solved up until this point either: 
                1) It is the first codon being used in this loop or 2) It differs from the previous ones in every position (or the 1st and last positions)"""
                count += 1                              
                state1 = 'S' + str(count)                       #new state name
                Delta.update({(origin, nucleotide1): state1})   #path from q0 new state with current nucleotide
                Q.add(state1)                                   #must update Q
                newstates1.append(state1)

                count += 1
                state2 = 'S' + str(count)                       #new state name
                Delta.update({(state1, nucleotide2): state2})   #path from previous state to new state with second nucleotide
                Q.add(state2)                                   #must update Q
                newstates2.append(state2)

                Delta.update({(state2, nucleotide3): destination})   #path from previous state to midpoint
        
            paths.append(i)

        for nucleotide in Sigma:
            if not Delta.has_key((origin,nucleotide)):
                Delta.update({(origin,nucleotide): sideA})

        for j in newstates1:
            for nucleotide in Sigma:
                if not Delta.has_key((j,nucleotide)):
                    Delta.update({(j,nucleotide): sideB})

        for j in newstates2:
            for nucleotide in Sigma:
                if not Delta.has_key((j,nucleotide)):
                    Delta.update({(j,nucleotide): origin})
        if n==2:
            for a in newstates1: F.add(a)
            for b in newstates2: F.add(b)


    for i in Sigma:
        Delta.update({('S0A', i): 'S0B'})
        Delta.update({('S0B', i): 'S0'})
        Delta.update({('SmidA', i): 'SmidB'})
        Delta.update({('SmidB', i): 'Smid'})
        Delta.update({('SendA', i): 'SendB'})
        Delta.update({('SendB', i): 'Send'}) 

    dfa = mk_dfa(Q,Sigma,Delta,q0,F)
    return dfa

"""Return the difference in nucleotides between the two codons - WORKS
    Return in the form: list(difference, n1, n2, n3) - if they differ in the first nucleotide n1 = 1. If they do not, n1=0, and so forth """
def get_difference(codon1, codon2):
    count = n1 = n2 = n3 = 0        
    if codon1[0] is not codon2[0]:  #compare the first nucleotide
        n1=1                        #if they are different, n1=1
        count += 1
    if codon1[1] is not codon2[1]:  #compare the second nucleotide
        n2=1                        #if they are different, n2=1
        count += 1
    if codon1[2] is not codon2[2]:  #compare the third nucleotide
        n3=1                        #if they are different, n3=1
        count += 1
    if count==2 and n1 == 1 and n3 == 1:  #if they differ in the first and last one, but not the second, they might as well differ in all nucleotides - for DFA purposes
        count = 3
    result = [count,n1,n2,n3]
    return result

""" generate_random_dna():  DOES NOT WORK YET!!!
Creates a random sequence of DNA
Type: Organism (Standard, Vertebrate Mitochondrial...)
Size: The number of nucleotides
"""
def random_dna(type, size):
    assert size >= 9, "There must be at least 9 nucleotides"
    assert size % 3 == 0, "Only codons!"
    ncodons = size / 3

    table = CodonTable.unambiguous_dna_by_id[type]
    start = table.start_codons
    stop = table.stop_codons
    dna = ''
    T1 = ''
    T2 = ''

    T1length = 1
    T2length = 1

    """S -> TCT means: S -> (T1)(C)(T2)"""
    
    C = rand_coding_seq(start,stop,(randint(3,ncodons))*3)
    Ccodons = len(C)/3

    T1length = randint(0,ncodons-Ccodons)
    while T1length > 0:
        curlength = randint(0,T1length)*3
        if randrange(0,2,1) == 0 and curlength >= 9: current = rand_coding_seq(start,stop,curlength)
        else: current = rand_noncoding_seq(start,curlength,True)
        T1length -= curlength/3
        T1 += current
    T1codons = len(T1)/3

    T2length = ncodons - Ccodons - T1codons
    while T2length > 0:
        curlength = randint(0,T2length)*3
        if randrange(0,2,1) == 0 and curlength >= 9: current = rand_coding_seq(start,stop,curlength)
        else: current = rand_noncoding_seq(start,curlength,True)
        T2length -= curlength/3
        T2 += current

    return Seq(str(T1 + C + T2), IUPAC.unambiguous_dna)

""" generate_coding_sequence(): WORKS?
Creates a random coding sequence of DNA
startcodons: do not allow the sequence to contain one of these
size: the number of nucleotides in the return sequence
onlycodons: whether or not the strand should be a multiple of three
"""
def rand_noncoding_seq(startcodons, size, onlycodons=False):
    if onlycodons:  assert size % 3 is 0, "Size must be a multiple of 3"
    assert isinstance(startcodons,list), "startcodons must be a list"
    assert len(startcodons)>0, "There must be at least one start codon"
    if size <= 0: return ''

    nucleotides = ['A','T','G','C'] #all possible nucleotides
    ncodons = size / 3
    remainder = size % 3
    strand = ''

    for n in range(0, ncodons):

        codon = startcodons[0]
        while startcodons.__contains__(codon):
            n1 = nucleotides[randint(0,3)]
            n2 = nucleotides[randint(0,3)]
            n3 = nucleotides[randint(0,3)]
            codon = n1 + n2 + n3    #create a codon from the three random nucleotides
        strand += codon #add the new codon to the strand

    """add additional nucleotides if necessary:"""
    if remainder == 1:
        strand += nucleotides[randint(0,3)]
    elif remainder == 2:
        strand += nucleotides[randint(0,3)]
        strand += nucleotides[randint(0,3)]

    return strand

""" generate_coding_sequence(): WORKS?
Creates a random coding sequence of DNA
startcodons: must start with one of these codons
stopcodons: must end with one of these codons
size: the number of nucleotides in the return sequence
"""
def rand_coding_seq(startcodons, stopcodons, size):
    assert size % 3 is 0, "Size must be a multiple of 3" #otherwise invalid coding sequence
    assert isinstance(startcodons,list) and isinstance(stopcodons,list), "startcodons and stop codons must both be lists"
    assert len(startcodons)>0 and len(stopcodons)>0, "There must be at least one start and one stop codon"
    if size <= 8: return ''

    nucleotides = ['A','T','G','C']
    ncodons = (size / 3) - 2  #must have a start codon at the beginning and stop at the end (thus minus 2 codons)

    strand = str(startcodons[randint(0,len(startcodons)-1)]) #begin strand with a start codon

    for n in range(0, ncodons):
        codon = stopcodons[0]
        while stopcodons.__contains__(codon):
            n1 = nucleotides[randint(0,3)]
            n2 = nucleotides[randint(0,3)]
            n3 = nucleotides[randint(0,3)]
            codon = n1 + n2 + n3    #create a codon from the three random nucleotides
        strand += codon #add the new codon to the strand

    strand += str(stopcodons[randint(0,len(stopcodons)-1)]) #end strand with a stop codon

    return strand

def validate_dna(strand, type): #Works
    dfa = mk_dna_dfa(type);
    return accepts(dfa,dfa['q0'],strand)

def test_functions(low, high):
    species = [1,2,3,4,5,6,9,10,11,12,13,14,15,16,21,22,23]
    for n in range(low, high):
        for m in range(1,17):
            type = species[m]
            dna = random_dna(type,n%3)
            assert len(dna)==n
            assert validate_dna(dna,type)==True
    return True


def mk_dfa(Q, Sigma, Delta, q0, F):

    assert len(Q)>0, "Q must not be empty"
    assert len(Sigma)>0, "Sigma must not be empty"

    for i in Sigma:
        assert isinstance(i, str), "Sigma must be a set of single character strings"
        assert len(i)==1, "Sigma must be a set of single character strings"

    assert Q.__contains__(q0), "q0 must belong to Q"

    assert isinstance(F, set), "F must be a non-empty set"
    assert len(F&Q)>0 and len(F-Q)==0, "F must be a subset of Q"

    assert isinstance(Delta, dict), "Delta is not a total function represented as a hash table"

    for i in Delta: 
        assert Q.__contains__(i[0]), "For key pairs (q,c) every q must be in Q"
        assert Sigma.__contains__(i[1]), "For key pairs (q,c) every c must be in Sigma"
        assert Q.__contains__(Delta[i]), "All values must be in Q"
        
    dfa = dict({"Q":Q, "Sigma":Sigma, "Delta":Delta, "q0":q0, "F":F})
    return dfa