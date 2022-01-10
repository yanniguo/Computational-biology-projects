import math
import random
from Bio import Entrez, SeqIO

""""
Title: Computational Biology 309 Lab 2 -  DosR binding site motif search in Mycobacterium tuberculosis

Authors: Xinduo Fan, Yanni Guo and Kairuo Yan

Date: 10/6/2021
"""

def getCounts(motifs):
    """Return nt counts given a list of motif strings."""
    
    counts = []
    for i in range(len(motifs[0])):
        counts.append({'A': 0, 'C': 0, 'G': 0, 'T': 0})
    
    for motif in motifs:
        motif = motif.upper()
        for i in range(len(motif)):
            counts[i][motif[i]] += 1
    
    return counts
    
def getProfile(motifs):
    """Get a profile from a set of motifs."""
    
    counts = getCounts(motifs)
    total = len(motifs)
    for i in range(len(counts)):
        for nt in 'ACGT':
            counts[i][nt] /= total
    return counts  # really now a profile
    
def getProfileLaplace(motifs):
    """Get a profile from a set of motifs, adding pseudocounts 
       (Laplace's rule of succession) to prevent zero-probability events.
    """
    
    counts = getCounts(motifs)
    total = len(motifs) + 4
    for i in range(len(counts)):
        for nt in 'ACGT':
            counts[i][nt] += 1
            counts[i][nt] /= total
    return counts  # really now a profile
    
def getMostProbable(sequence, profile, k):
    """Return the most probable k-mer in sequence, given a profile."""
    
    maxProb = -1
    best_kmer = ''
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i + k]
        prob = 1
        for j in range(k):
            prob *= profile[j][kmer[j]]
        if prob > maxProb:
            maxProb = prob
            best_kmer = kmer
    return best_kmer
    
def getMostProbableMotifs(sequences, profile, k):
    """Return the most probable set of k-mers from a set of sequences,
       given a profile."""
       
    motifs = []
    for sequence in sequences:
        motifs.append(getMostProbable(sequence, profile, k))
    return motifs
    
def getConsensus(motifs):
    """Return the consensus sequence for a set of motifs."""
    
    counts = getCounts(motifs)
    
    consensus = ''
    for i in range(len(counts)):
        majority = 'A'
        for nt in counts[i]:
            if counts[i][nt] > counts[i][majority]:
                majority = nt
        consensus += majority
        
    return consensus
    
def getScore(motifs):
    """Return the score for a set of motifs."""
    
    counts = getCounts(motifs)
    consensus = getConsensus(motifs)
    
    t = len(motifs)
    score = 0
    for i in range(len(consensus)):
        nt = consensus[i]
        diff = t - counts[i][nt]
        score += diff
        
    return score



def readSequence(name):
    """
    read the sequences in a fasta file and return the list with each sequence as a string value
    
    Parameters
    ----------
    name : (string) the name of the FASTA file that contains the sequences

    Returns
    -------
    sequences: (list of strings)  list with each sequence as a string value 
    """

    record = SeqIO.index( (name + '.fasta'),'fasta')      #index function reads the dna sequence into a dictionary then convert to a list of strings
    sequences = []                                        #empty list to store sequence strings  
    for id in record:                                     #for loop iterates dna sequence and add sequence into a list of strings  
        sequences.append(str(record[id].seq))

    return sequences                                      #return list of sequences 
                          

def GibbsSampling(dna, k, t, N):
    """
    randomly choose a sequence and adjust our choice of motif from 
    it with Gibbs Sampling technique for N times
    
    Parameters:
    ------------
        dna: (list of strings) the list of sequences; 
        k: (int) the length of the motifs we’re looking for; 
        t: (int) the number of sequences which is same as the number of motifs we are expected to find; 
        N: (int) the number of times we randomly choose a sequence and adjust our choice of motif from it with Gibbs Sampling technique

    Returns:
    ----------
    BestMotifs:(list of strings) the list of motifs we found with the lowest score
    """
    
    motifs = []                                                             #empty list to store motifs
    for i in range(t):                                                      #iterate over dna sequence list
        index = random.randrange(len(dna[i])-k+1)                           #randomly choose a motif from each sequence to create the first set of motifs
        motifs.append(dna[i][index:index+k])
    BestMotifs = motifs                                                     #records the randomized set of motifs
    BestScore = getScore(motifs)                                            #records the score of the motifs
    for j in range(N):                                                      #repeat process to adjust motif choices N times
        i = random.randrange(t)                                             #randomly select a motif from the current set
        del motifs[i]                                                       #remove the selected motif from the set
        profile = getProfileLaplace(motifs)                                 #call function to create a profile from current motif set without the selected motif
        prossibilities = {}                                                 #empty dictionary to store possible motifs in the ith sequence
        sum = 0                                                             #tracks the sum of possibilities
        for x in range(len(dna[i])-k+1):                                    #iterate over dna sequence to calculate the possibility of each motif canditate (using current profile) and the sum of these prosibilities
            candidate = dna[i][x:x+k]                                       
            if candidate not in prossibilities:
                prossibilities[candidate] = 0
                for position in range(k):                                   
                    nt = candidate[position]
                    prossibilities[candidate] += profile[position][nt]
            sum += prossibilities[candidate]
        pointer = 0                                                         #used to iterate through the ith sequence
        possibility_keeper = 0                                              #tracks current possibility's position relative to the sum of all posibilities
        chosen = random.random()                                            #generate a random number lower than the sum
        chosen *= sum
        while (chosen > possibility_keeper):                                #go through each motif in the ith sequence to find where the random number (chosen variable) locates in the possibility map
            candidate = dna[i][pointer:pointer+k]
            possibility_keeper += prossibilities[candidate]
            pointer += 1                                                    #end at the motif which the random number (chosen variable) points to
        motifs.insert(i,candidate)                                          #insert the chosen motif back into the motif set at the ith position
        score = getScore(motifs)                                            #calculate score for current set of motifs
        if score < BestScore:                                               #update BestMotifs and BestScore to the new current set of motifs if the current score is lower than the BestScore
            BestMotifs = motifs
            BestScore = score
    motif_score = {}
    motifs = BestMotifs
    for i in range(t):
        motif = motifs[i]
        del motifs[i]
        if motif not in motif_score:
            motif_score[motif] = 0
        score = BestScore - getScore(motifs)
        motif_score[motif] = score
        motifs.insert(i,motif)
    while len(BestMotifs) > 25:
        motif = max(motif_score, key=motif_score.get)
        BestMotifs.remove(motif)
        del motif_score[motif]
    return BestMotifs                                                       #return best set of motifs and its score after N iterations are completed

def repeatGS(sequences,times, N):
    """
    repeatly call the GibbsSampling function, and keep track of the best set of motifs with the lowest score
    
    Parameters:
    ------------
        sequences: (list of strings) the list of sequences; 
        times : (int) the number of times we call the GibbsSampling function
        N: (int) the number of times we randomly choose a sequence and adjust our choice of 
                 motif from it with Gibbs Sampling technique in the GibbsSampling function
                 
    Returns:
    -----------
    BestM:(list of strings) the list of motifs we found with the lowest score
    """
    
    bestM = GibbsSampling(sequences, 20, len(sequences), N)      #initialize variable that calls GibbsSampling function for best motifs
    bestS = getScore(bestM)                                      #initialize variable that calls getScore function for score of best motifs
    for i in range(times-1):                                     #repeativley call GibbsSampling function and getScore function
        M = GibbsSampling(sequences, 20, len(sequences), N)
        S = getScore(M)
        if S<bestS :                                             #update bestM and bestS variables when a motif set gets a lower score than the current bestS            
            bestM = M 
            bestS = S
    return(bestM)                                                #return best set of motifs




def main():
    """
    Functions are called when given the upstream sequence of the MTB gene.
    Returns the best set of motifs, motifs scores, and consensus sequences
    """
    #section 2
    sequences = readSequence('upstream250')
   
    #section 1
    bestM = GibbsSampling(sequences, 20, len(sequences),1000)
    bestS = getScore(bestM)
    print(bestM, bestS)

    #section 3
    bestM = repeatGS(sequences, 1000,100)
    bestS = getScore(bestM)
    print(bestM, bestS)
    
    #section 4 
    Consensus = getConsensus(bestM)
    print(Consensus)

    #section 6
    motifs = ['TTCGTGACCGACGTCCCCAG','TTGGGGACTTCCGGCCCTAA','GCCGGGACTTCAGGCCCTAT','CATGGGACTTTCGGCCCTGT','GAGGGGACTTTTGGCCACCG','CCAGGGACCTAATTCCATAT','TTGAGGACCTTCGGCCCCAC','CTGGGGACCGAAGTCCCCGG','TTAGGGACCATCGCCTCCTG','TGGATGACTTACGGCCCTGA','TTGGGGACTAAAGCCTCATG','TCGGGGACTTCTGTCCCTAG','TTGGGGACCATTGACCCTGT','TTGAGGACCTAAGCCCGTTG','CACGGGTCAAACGACCCTAG','GGCGGGACGTAAGTCCCTAA','GAAGTGACGAAAGACCCCAG','CGGAGGACCTTTGGCCCTGC','GTGGGGACCAACGCCCCTGG']
    Consensus = getConsensus(motifs)
    score = getScore(motifs)
    print(Consensus,score)

main()