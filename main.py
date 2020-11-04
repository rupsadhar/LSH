import os
import pandas as pd
import binascii as bin
import numpy as np
import sympy as sympy
import random as random
import math
import time
data_list=[]
INF = 2**32
number_of_hash=100
threshold = 0.6
rows_of_B = 5
shingle_size = 7
docAsShingleSets={}

def parsing_data(inputQuery):
    """ Parses data form txt file and converts it to a list.

    :return: doc_list: list of documents list from the dataset
    """
    data = pd.read_csv("./dna_data/chimp_data-noN.txt", sep = "	")

    #each sequence stored as a list from the corpus data
    doc_list = data["sequence"].tolist()
    doc_list.append(inputQuery)
    return doc_list

def build_kmers(sequence, ksize):
    """ Builds k-shingles for each sequence.
    
    :param sequence: sequence of DNA per document
    :param ksize: Size of shingles
    :return: list of shingles of size k from each document 
    """
    kmers = []
    n_kmers = len(sequence) - ksize + 1

    for i in range(n_kmers):
        kmer = sequence[i:i + ksize]
        kmers.append(kmer)

    return kmers


shingles={}
docIDlist=set()

def shingle(inputQuery):
    """ Function to shingles documnets and find posting lists.

    :param inputQuery: The query input from the user
    :return: Document-wise Shingles, All unique shingles from corpus, The Posting Lists per shingle
    :rtype: Dictionary, set, Dictionary of sets
    """
    sequences = parsing_data(inputQuery)
    cnt=0
    shinglesInDocWords = set()
    PostingDict = {}
    print("\n Shingling the documents ")
    for i in sequences:

       
        docIDlist.add(cnt+1)
        shingle_list=build_kmers(i,shingle_size)

        # docAsShingleSets[cnt]=shingle_list
        templist=[]
        l=0
        for shingle_word in shingle_list:
            hashed_shingle = bin.crc32(shingle_word.encode('ASCII')) & 0xffffffff
            if hashed_shingle not in shinglesInDocWords:
                shinglesInDocWords.add(hashed_shingle)
            if hashed_shingle not in templist:
                templist.append(hashed_shingle)
                l+=1

            # Creating the Posting list for each shingle
            if hashed_shingle not in PostingDict:
                PostingDict[hashed_shingle] = set()  

            PostingDict[hashed_shingle].add(cnt)

        docAsShingleSets[cnt]=templist
        cnt+=1
    #print(cnt)

    list_of_unique_shingles=[]
    count=0
    for i in shinglesInDocWords:
        list_of_unique_shingles.insert(count,i)
        count+=1


    return docAsShingleSets, list_of_unique_shingles, PostingDict, docIDlist


def matrixGenerator(allShingles, invertedIndexTable):
    """ This function indexes the shingles and returns a boolean matrix of shingle versus document.

    :param allShingles: list of all shingles in the corpus
    :param invertedIndexTable: dictionary with posting list for each shingle
    :return: Matrix 
    :rtype: 2D Array
    """
    posting_list_1=[]
    matrix_with_index = {}
    index = 0
    # indexing the shingles
    print("\n Generating Boolean Matrix\n")
    for shingle in allShingles:
        matrix_with_index[shingle] = index
        index += 1

    # shingle document matrix
    check=0
    matrix = np.zeros([len(allShingles), 1681], dtype=int)
    for shingle in allShingles:
        posting_list_1 = invertedIndexTable[shingle]
        for d in posting_list_1:
            matrix[matrix_with_index[shingle]][int(d)] = 1  # Boolean value true for that document corresponding to a shingle
            check+=1
    return matrix

def pickRandomCoeffs(k, maxval):
    ''' This function returns k number of unique random values.

    :param k: Number of random values
    :param maxval: maximum value for randint
    :return: A list of random integers
    '''
  # Create a list of 'k' random values.
    rand_no_list = []

    while k > 0:
    # Get a random shingle ID.
        randIndex = random.randint(0, maxval)

    # Ensure that each random number is unique.
        while randIndex in rand_no_list:
            randIndex = random.randint(0, maxval)

        rand_no_list.append(randIndex)
        k = k - 1

    return rand_no_list

def find_sign_matrix(matrix, numOfShingles):
    """ This function picks two random coefficient values and genrates a hashfunction of form h(x)=(ax+b)%c
    All values are initialised to infinities and each row is mapped to lowest hash function that has a
    boolean true for that shingle. This new matrix called sigmatrix is returned.
    example
    matrix= 
    [[1, 0, 0, 1],
    [0, 0, 1, 0],
    [0, 1, 0, 1],
    [1, 0, 1, 1],
    [0, 0, 1, 0]]
    
    coeffA = [1,1]
    coeffB = [1,3]
    c = 5
    required output is [[1, 3, 0, 1], [0, 2, 0, 0]].

    :param matrix: boolean matrix of shingles vs docs
    :param numOfShingles: total number of shingles in corpus
    :return: The signature matrix based on the hash values
    """
    print("\n Generating signature Matrix")
    c = numOfShingles
    while not sympy.isprime(c):
        c += 1

    coefficient_1 = pickRandomCoeffs(number_of_hash, numOfShingles-1)
    coefficient_2 = pickRandomCoeffs(number_of_hash, numOfShingles-1)

    rows, cols, sigrows = len(matrix), len(matrix[0]), len(coefficient_1)
    # initialize signature matrix with maxint
    signature_Mat = np.full((sigrows, cols), INF)
    

    for r in range(rows):
        hashvalue = []
        for h in range(sigrows):
            hashvalue.append((coefficient_1[h] + coefficient_2[h]*r) % c)# Hash each row
        # if data != 0 and signature > hash value, replace signature with hash value
        for col in range(cols):
            if matrix[r][col] == 0:
                continue
            for i in range(sigrows):
                if signature_Mat[i, col] > hashvalue[i]:
                    signature_Mat[i, col] = hashvalue[i]
    print("Signature matrix")
    print(signature_Mat)
    return signature_Mat

def getbestb(threshold,number_of_hash, eps=1e0):
    """ Calculates the best value for no. of bands.

    :param threshold: difined threshold
    :param number_of_hash: number of hash functions
    :param eps: Constant value
    :return: The best value for b by solving an equation
    """
    for b in range(1, number_of_hash+1):
        opt = b*math.log10(b)
        val = -1 * number_of_hash * math.log10(threshold)
        if opt > val-eps and opt < val+eps:
            print("\n Processing LSH with number of bands : %d" % (np.round(b)))
            return np.round(b)


def lsh(B_BANDS, docIdList, sig):
    """ Applies the LSH algorithm. This function first divides the signature matrix into bands and hashes each column onto buckets.

    :param B_BANDS: Number of bands in signature matrix
    :param docIdList: List of document ids
    :param sig: signature matrix
    :return: List of document to its hash along with the buckets
    """
    numHash = number_of_hash
    bands = getbestb(threshold,number_of_hash)
    rows = numHash / bands

    d = 1681
    # Array of dictionaries, each dictionary is for each band which will hold buckets for hashed vectors in that band
    buckets = np.full(bands, {})
    # Mapping from docid to h to find the buckets in which document with docid was hashed
    docth = np.zeros((d, bands), dtype=int)  # doc to hash
    for i in range(bands):
        for j in range(d):
            low = int(i*rows) # First row in a band
            high = min(int((i+1)*rows), numHash)# Last row in current band
            l = []
            for x in range(low, high):
                l.append(sig[x, j])  # Append each row into l
            h = int(hash(tuple(l))) % (d+1)
            try:
                buckets[i][h].append(j) # If a bucket corresponds to this hash value append this document into it
            except:
                buckets[i][h] = {j}
            docth[j][i] = h
    # print(docth)
    return docth, buckets
def jacsim(doc1, doc2, docsAsShingleSets,sign_matrix):
    """ Jaccard Similarity.

    :param doc1: First doc to be compared
    :param doc2: Second doc to be comapred
    :param docsAsShingleSets: Document wise shingles
    :param sign_matrix: The Signature matrix
    :return: The jaccard similarity value
    """
    document1 = sign_matrix[:,doc1]
    document2 = sign_matrix[:,doc2]
    intersect = sum(bool(x) for x in np.logical_and(document1, document2))
    return (intersect / len(document1))

def get_similar(dn,docIdList,buckets,docth,docsAsShingleSets,sign_matrix):
    """" This function finds similar documents given a query document after hashing and bucketing the query document
    It also evaluates based on various similarity criterion, namely, Jacard similarity, Euclidean distance
    and cosine similarity.

    :param dn: The query document number
    :param docIdList: List of doc ids
    :param buckets: List of buckets
    :param docth: doc to hash list
    :param docAsShingleSets: Document wise shingles
    :param sign_matrix: The Signature matrix
    :return: List of similar documents based on decreasing similarity amount
    """
    if dn not in docIdList:
        raise KeyError('No document with the given name found in the corpus.')

    docid = dn
    # Collection of documents similar to docid
    coooo = []
    # taking union of all buckets in which docid is present
    for b, h in enumerate(docth[docid]):
        coooo.extend(buckets[b][h])
    coooo = set(coooo)
    print("\nComparing with docs")
    print(coooo)

    # Similar documents
    sim_list = []
    for doc in coooo:
        if doc == docid:
            continue
        sim = jacsim(docid, doc, docsAsShingleSets,sign_matrix)
        sim_list.append((sim, doc))
    sim_list.sort(reverse=True)
    return sim_list


def main():
    """ Main Function.
    """
    inputQuery = input("Enter query string:")
    tic = time.time()
    docsAsShingleSets, allShingles, PostingDict, docIDlist = shingle(inputQuery)
    toc = time.time()
    print("Time taken = ", toc-tic)

    tic = time.time()
    matrix = matrixGenerator(allShingles,PostingDict)
    print(matrix)
    toc = time.time()
    print("Time taken = ", toc-tic)

    tic = time.time()
    sign_matrix = find_sign_matrix(matrix,len(allShingles))
    toc = time.time()
    print("Time taken = ", toc-tic)

    tic = time.time()
    BANDS=20
    docth,buckets = lsh(BANDS,docIDlist,sign_matrix)
    toc = time.time()
    print("Time taken = ", toc-tic)

    query_id = len(docAsShingleSets)-1
    inputDocID=query_id

    sim_docs = get_similar(int(inputDocID),docIDlist,buckets,docth,docsAsShingleSets,sign_matrix)

    print("Calculating Jaccard similarities....\n")


    found = 0
    for sim, doc in sim_docs:
        if sim >= threshold:
            found = 1
            print('Document Name: ' + str(doc), 'Similarity: ' + str(sim) + '\n')
    if found == 0:
        print("NO similar docs for the given threshold")

if __name__ == "__main__":
    main()