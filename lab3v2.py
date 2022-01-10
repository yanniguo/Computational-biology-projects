import math
import random
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import queue

def readpairs(name,k):
    """
    Read the FASTA file and return a list of reads and a list of k-mers from the reads
    
    Parameters
    ----------
    name : FASTA file of bacterial genome
    k : (integer) length of reads

    Returns
    -------
    A list of reads and a list of k-mers from the read

    """
    kmers = []                                                  #empty list to store k-mers
    sequences = []                                              #empty list to stores reads
    for record in SeqIO.parse( (name + '.fasta'),'fasta'):      #iterate over FASTA file and parse function separates read labels and read-pairs
        sequences.append(str(record.seq))                       #add reads to reads list
    for i in range( len(sequences) ):                           #iterate over index of reads list to generate k-mers
        for j in range( len(sequences[i])-k+1 ) :
            kmer = sequences[i][j:j+k]                          
            kmers.append(kmer)                                  #add k-mers to k-mers liist
    return (sequences,kmers)                                    #return lists for reads and k-mers


def create_deBruijn(kmers):
    """
    Create a de Bruijn graph from a set of k-mers
    
    Parameters
    ----------
    kmers : (string) pattern of bases from the readsof the bacteria genome

    Returns
    -------
    Return a graph dictionary, a dictionary for in-degrees of each node, and a 
    dictionary for out-degress of each node

    """
    count = {}                                    #empty dictionary to record k-mers
    graph = {}                                    #empty dictionary to record nodes and list of nodes connected by a directed edge
    in_degrees = {}                               #empty dictionary to record in-degrees of nodes
    out_degrees = {}                              #empty dictionary to record out-degrees of nodes
    for i in range(len(kmers)):                   #for loop iterates over list of k-mers  
        if kmers[i] not in count:                 #if k-mers index not in k-mers dictionary, then node prefix is in the k-mer backward one position  
            node_p = kmers[i][:-1]
            if node_p not in graph:               #if node prefix not in graph dictionary then node prefix is in empty list      
                graph[node_p] = []
                in_degrees[node_p] = 0            #initialize in-degree of node prefix 
                out_degrees[node_p] = 0           #initialize out-degree of node prefix 
            node_s = kmers[i][1:]                 #node suffix equals the k-mers index position
            if node_s not in graph:               #if node suffix not in graph dictionary then node suffix is in empty list  
                graph[node_s] = []
                in_degrees[node_s] = 0            #initialize in-degree of node suffix
                out_degrees[node_s] = 0           #initialize out-degree of node suffix 
            graph[node_p].append(node_s)          #add node prefix to node suffix list  
            in_degrees[node_s] += 1               #in_degree of node suffix add one value  
            out_degrees[node_p] += 1              #out-degree of node prefix add one value
            count[kmers[i]] = 0                   #k-mers index in dictionary equals 0  
    return (graph,in_degrees,out_degrees)         #return graph, in-degree, and our-degree dictionaries   


def GetMaximalNonbranchingPaths(graph, inDegrees, outDegrees):
    """
    Find all maximal non-branching paths in the de Bruijn graph. Which are paths with internal
    nodes that have in-degree and out-degree equal to one except for the starting and ending node
    
    Parameters
    ----------
    graph : dictionary with keys as the nodes ((k-1)mers) and values as a list of directed edge
            that connects nodes and keys
    inDegrees : dictionary of keys as nodes and values as the corresponding number of in-degrees
    outDegrees : dictionary of keys as nodes and values as the corresponding number of out-degrees
        
    Returns
    -------
    paths: list of maximal non-branching paths in the de Bruijn graph, consist of sets of k-mers

    """
    paths = []                                              #empty list to store maximal non-branching paths
    nodes11 = {}                                            #empty dictionary to store set of nodes with in-degree and out-degree both equal to 1
    for id in graph:                                        #for loop iterate over index of graph
        if inDegrees[id] == 1 and outDegrees[id] == 1:      #if in-degree and out-degree of index in graph equals 1
            nodes11[id] = 0                                 #then the values of nodes11 dictionary equals 0
    
    for start in graph:                                     #iterate over each node start in graph    
        if start not in nodes11:                            #if node start is not in nodes11 dictionary
            for current in graph[start]:                    #iterate over each current node in adjacent list of start
                path = [start,current]                      #path equals node start and node current
                while (current in nodes11):                 #while current node in is nodes11 dictionary
                    next_ = graph[current][0]               #set variable for node next to current node
                    path.append(next_)                      #append next node to paths list
                    graph[current].remove(next_)            #remove next node from current node list
                    current = next_                         #set current node as next node
                paths.append(path)                          #append paths to list of paths
    
    for start in nodes11:                                   #iterate over each start node in the nodes11 dictionary
        if (graph[start] != []):                            #if start node does not have an adjacent node then follow conditions
            path = [start]                                  #set path equal to start node
            current = graph[start][0]                       #set current node to the node adjacent to start node
            while (current != start):                       #while current node doesn't equal start node add current node to path
                path.append(current)                    
                next_ = graph[current][0]                   #set next variable equal to the node adjacent to the current node
                graph[current].remove(next_)                #remove current and next node from the graph
                current = next_                             #set current node equal to the next node
            paths.append(path)                              #add path to paths list
    
    return paths                                            #return paths list

def constructContigs(paths):
    """
    Function creates contigs (overlapping DNA bases) from a set of paths

    Parameters
    ----------
    paths : list of kmers that is the maximal non-branching path in the de Bruijn graph

    Returns
    -------
    Returns a list of contigs

    """
    contigs = []                                #empty list to store contigs
    for i in range(len(paths)):                 #iterate over index of paths list
        path = paths[i]                         #set variable path equal to index of paths list
        contig = path[0]                        #set variable contig equal to first index of path
        for j in range(1,len(path)):            #iterate over next index of path
            contig = contig + path[j][-1]       #set contig equal to sum of contig and path index minus 1
        contigs.append(contig)                  #add contig to contigs list
    return (contigs)                            #return contigs as a list 

def calculateN50(l):
    """
    Function used N50 statistic to measure the quality of the assembly of contigs
    
    Parameters
    ----------
    l = (list) set of reads
    
    Returns
    ---------
    Returns N50 value: maximal contig length that the contigs greater than or equal to the length
    that is equal to half of the sun of lengths of contigs
    """
    
    middle = sum(l)//2     #set middle variable to calculate half of sum of the lengths of all contigs
    l.sort()               #sort the list of reads in order from minimum to maximum  
    i = -1                 #set index of congs equal to -1         
    N50 = l[i]             #set N50 variable equal to reads list
    count = N50            #set count varaible equal to N50
    while count < middle:  #while count is less than 50% of the total lengths of contigs
        i -= 1             #index equals minus 1 value
        N50 = l[i]         #set N50 equal to index of reads list
        count += N50       #set count equal to add value of N50
    return N50             #return N50 value

def orderingContigs(contigs, reads):
    """
    Function determine ordering of contigs with pair-end reads of original FASTA file

    Parameters
    ----------
    contigs : list of contigs
    reads : (string) reads of origial FASTA file

    Returns
    -------
    Returns new FASTA file called "lab3_result.fasta"

    """
    order_graph = {}                                                    #empty dictionary to record ordering of 
    indegree = {}
    final_contigs = []
    for i in range( len(contigs) ):                                     #iterate over the length of contigs
        current = contigs[i]                                            #set the ith contig as current
        for j in range( 0, len(reads), 2 ):                             #the second index of the read label
            if reads[j] in current:
                if current not in order_graph:                          #if ith contig is not in order_graph dictionary
                    order_graph[current] = []                           #set the ith contig as key and assign it with an empty list
                if current not in indegree:                             #if ith contig is not in indegree dictionary
                    indegree[current] = 0                               #set the indgree of the contig as 0
                for x in range( len(contigs) ):
                    if i != x:
                        if reads[j+1] in contigs[x]:                     #if the reads is not in contigs
                            follow = contigs[x]                          #set the contig as follow
                            if follow not in indegree:                   #set the indegree of follow as 0   
                                indegree[follow] = 0
                            order_graph[current].append(contigs[x])      #combine contigs
                            indegree[follow] += 1                        #add indegree 1 to the follow variable
    
    #order = queue.Queue()
    order = []
    for key in order_graph:                         #iterate over the order_graph dictionary we set before
        if indegree[key] == 0:                      #if the indegree of the contig is 0
            #order.put(id)
            order.append(key)                       #append the first contig that with the indegree of 1 to order list

    while (len(order)!=0):  #not(order.empty())
        key = order.pop(-1)
        final_contigs.append(key)
        if key in order_graph:
            succeessor = order_graph[key]
        if (len(succeessor) != 0 ):                #if the indegree of the key(contig) is not 0  
            for i in range( len(succeessor) ):
                node = succeessor[i]               #node is the contig
                indegree[node] -= 1                #the indegree of node minus 1
                if (indegree[node] == 0):
                    order.append(node)             #append all nodes that with indegree of 0 into the order 
    #print(len(final_contigs))
    fasta_seq = []
    for i in range(len(final_contigs)):
        contig = final_contigs[i]                                        #set the contig as each key in final list
        record = SeqRecord(Seq(contig), "contig%i" % (i + 1), "", "")    #build the fasta file for contigs  
        fasta_seq.append(record)                                         #append the record to fasta file
    j = 1
    for i in range(len(contigs)):
        if contigs[i] not in final_contigs:
            contig = contigs[i]                                                             #set the contig as each key in contigs liset
            record = SeqRecord(Seq(contig), "contig%i" % ( len(final_contigs) +j ), "", "") #build the fasta file for contigs
            fasta_seq.append(record)                                                        #append the record to fasta file 
            j+=1 
    SeqIO.write(fasta_seq, "lab3_result.fasta", "fasta")                                    #write out FASTA file


def main():
    """
    Functions are called when given the FASTA file for C_ruddii bacterium.
    Returns a list of reads and a list of k-mers, graph dictionary, in-degree dictionary, out-degree dictionary,
    list of contigs, best N50, best K value, number of contigsm, and the length of each contig.

    """
    num_of_contigs = []         #empty list to store number of contigs
    contig_max_len = []         #empty list to store maximum length of contig
    Best_N50 = 0                
    Best_K = 20
    Best_set_contig_len = []    #empty list to store best set of contig length
    num_of_contigs = 0
    Best_contigs = []
    for j in range(20,41):                        #between 20 and 41 as values for k, variables are assigned to call functions
        reads,kmers = readpairs('C_ruddii',j)
        graph,in_degrees,out_degrees = create_deBruijn(kmers)
        paths = GetMaximalNonbranchingPaths(graph, in_degrees, out_degrees)
        contigs = constructContigs(paths)
        a = []                                   #empty list to store length of contigs       
        for i in range( len(contigs) ):          #iterate over length of contigs and add them to list a
            a.append(len(contigs[i]) )
        N50 = calculateN50(a)                    #calculate N50 for maximal contig length 
        if N50 > Best_N50:
            Best_N50 = N50
            Best_K = 20+i
            num_of_contigs = len(contigs)
            Best_set_contig_len = a
            Best_contigs = contigs
            
    print('Highest N50 value:',Best_N50,'when K=', Best_K,'the set of contigs has',num_of_contigs, 'contigs each has a length of', Best_set_contig_len)
    #print(paths)
    orderingContigs(Best_contigs,reads)
    
main()

"""
print(num_of_contigs,contig_max_len)
print(min(num_of_contigs),num_of_contigs.index(min(num_of_contigs)))
print(max(contig_max_len), contig_max_len.index(max(contig_max_len)))
"""