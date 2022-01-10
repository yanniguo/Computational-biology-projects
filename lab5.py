from Bio import Entrez, SeqIO, Phylo
from io import StringIO
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
#import heatmap

def readSequence(name,sequences):
    """
    read the sequences in a fasta file and return the list with each sequence as a string value
    
    Parameters
    ----------
    name : (string) the name of the FASTA file that contains the sequences
    sequences: (list of strings)  list with each sequence as a string value

    Returns
    -------
    Returns list of sequences as strings
    """

    record = SeqIO.index( (name + '.fasta'),'fasta')      #index function reads the dna sequence 
    for id in record:                                     #for loop iterates dna sequence and add sequence into the given dictionary
        if id not in sequences:
            sequences[id] = str(record[id].seq)

    return sequences                              

def propotion_B117_variants(sequences,labels,is_1205 = False):
    
    labels_for_theday = []
    N501_pos = 23063
    H69_70_pos = 21765
    P681_pos = 23604
    if is_1205:
        N501_pos = 23071
        H69_70_pos = 21764
        P681_pos = 23612
    count = 0
    for key in sequences:
        if (sequences[key][N501_pos-1] == 'T') and (sequences[key][H69_70_pos-1:H69_70_pos+5] == '------') and (sequences[key][P681_pos-1] == 'A'):
            count+=1
            labels_for_theday.append(key)
    labels.append(labels_for_theday)
    propotion = count / ( len(sequences) - 1 )
    return (propotion, labels )

def excercise_7():
    name = ['11-03','11-10','11-27','12-05','12-08']
    sequences = {}
    labels = []
    propotions = []

    sequences = readSequence(name[0]+' alignments',sequences)
    propotion_11_03 , labels= propotion_B117_variants(sequences,labels)
    propotions.append(propotion_11_03)

    sequences = readSequence(name[1]+' alignments',sequences)
    propotion_11_10 , labels= propotion_B117_variants(sequences,labels)
    propotions.append(propotion_11_10)

    sequences = readSequence(name[2]+' alignments',sequences)
    propotion_11_27 , labels= propotion_B117_variants(sequences,labels)
    propotions.append(propotion_11_27)

    sequences = readSequence(name[3]+' alignments',sequences)
    propotion_12_05 , labels= propotion_B117_variants(sequences,labels,True)
    propotions.append(propotion_12_05)

    sequences = readSequence(name[4]+' alignments',sequences)
    propotion_12_08 , labels= propotion_B117_variants(sequences,labels)
    propotions.append(propotion_12_08)



    plt.figure()
    ax = plt.bar(name,propotions)
    plt.bar_label(ax)
    plt.title('propotion of genomes with B117 variants for each date')
    plt.show()

    return(propotions,labels)

def get_file_lines(filepath):
    lineStrings = []
    FILE = open(filepath, "r")
    for line in FILE:
        if (line[-1] == "\n"): lineStrings.append(line[:-1])
        else: lineStrings.append(line)
    FILE.close()
    return lineStrings

def lines_to_dataframe(lines):
    names = []
    dist_dict = dict()
    for line in lines:

        linarr = line.split()
        percents = list(map(lambda x : float(x), linarr[2:]))
        anno = linarr[1].split("|")[1]

        if anno in dist_dict.keys(): print(anno)
        dist_dict[anno] = percents
        names.append(anno)

    DF = pd.DataFrame(columns = names)
    for sample in names: DF.loc[sample] = dist_dict[sample]

    return DF

def pim_to_dataframe(filepath : str) -> pd.DataFrame:
    lines = get_file_lines(filepath)
    return lines_to_dataframe(lines)

def proprocess_D(D):
    for i in range(len(D)):
        for j in range(len(D)):
            D.iloc[i,j] = 100 - D.iloc[i,j]
    return D
            
def computeM(D,taxa):
    n = len(D)
    R = np.zeros(n)                  # create a 1-D array of n zeros
    M = np.zeros((n, n))             # create an n x n array of zeros; refer to an element with M[i, j]

    for i in range(len(taxa)):
        for j in range(len(taxa)):
            R[i] += D.iloc[i,j]
    
    for i in range(n):
        for j in range(n):
            if i != j:
                M[i,j] = (n-2) * D.iloc[i,j] - R[i] - R[j]
    return M,R

def update_D(D,R,index,taxa):
    a = index[0]
    b = index[1]
    a_name = taxa[a]
    b_name = taxa[b]
    D_v = np.zeros(len(D)+1)
    outliners = []
    for k in range(len(taxa)):
        if (k!=a) and (k!=b):
            D_vk = (D.iloc[a,k] + D.iloc[b,k] - D.iloc[a,b]) /2
            D_v[k] = D_vk
    
    L_Va = ( (D.iloc[a,b] / 2) + ((R[a]-R[b])/ (2* (len(D)-2) )) )
    L_Vb = ( (D.iloc[a,b] / 2) + ((R[b]-R[a])/ (2* (len(D)-2) )) )

    if (L_Va>1):
        outliners.append(a_name)
    if (L_Vb>1):
        outliners.append(b_name)

    V_name = '('+ a_name +':'+ str(L_Va) +','+ b_name  + ':'+ str(L_Vb) +')' #'('+ a_name +','+ b_name +')'
    
    D.insert(len(D), V_name, D_v[:-1])
    D.loc[V_name] = D_v
    D = D.drop(columns = [taxa[a],taxa[b]])
    D = D.drop(index = [taxa[a],taxa[b]])
    taxa.remove(a_name)
    taxa.remove(b_name)
    taxa.append(V_name)

    return(D,taxa, V_name,outliners)

def read_martix(name):
    DF = pim_to_dataframe('./'+name)
    return (DF)

def nj(dist, taxa):
    dist = proprocess_D(dist)
    outliners = []
    while len(dist)>2:
        n = len(dist)
        M,R = computeM(dist,taxa)
        min_flat_index = np.argmin(M)        # get the "flat" index of the minimum value in array M
        (a, b) = np.unravel_index(min_flat_index, (n, n)) # convert a "flat" index into 2-D indices
        dist, taxa, newick, outliner = update_D(dist,R, (a,b), taxa)
        outliners.append(outliner)
    newick = '('+ taxa[0] +',' + taxa[1] + ':' + str(dist.iloc[0,1]) +')'
    return(newick,outliners)

def excercise_9():
    D_1103 = read_martix('HCOV19-ENGLAND-081220-D.pim')
    newick,outliner = nj(D_1103,list(D_1103.columns))
    T = Phylo.read(StringIO(newick), 'newick')
    Phylo.draw(T)
    return (outliner)
 
def main():
    propotions,labels = excercise_7()
    sample_id = []
    for label in labels[-1]:
        linarr = label.split()
        anno = linarr[0].split("|")[2]
        sample_id.append(anno)
    
    print(sample_id)
    outliner = excercise_9()
    outliners = ['EPI_ISL_816270', 'EPI_ISL_741303','EPI_ISL_733688','EPI_ISL_819570', 'EPI_ISL_736686', 'EPI_ISL_843006', 'EPI_ISL_741289','EPI_ISL_816861','EPI_ISL_702633','EPI_ISL_727871', 'EPI_ISL_816263','EPI_ISL_782764','EPI_ISL_813701','EPI_ISL_842663','EPI_ISL_733856']
    count = 0
    for i in outliners:
        for j in sample_id:
            if i == j:
                count += 1 
    print(count/len(outliners))


main()