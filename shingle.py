
import os
import json
import pandas as pd
import binascii as bin
import numpy as np
import sympy as sympy
import random as random
data_list=[]
INF = 2**32
NUM_HASH_FUNCS=150
docAsShingleSets={}
#parsing_data() returns the list of documents list from the dataset
def parsing_data() :

    data = pd.read_csv("./dna_data/chimp_data-noN.txt", sep = "	")

    #each sequence stored as a list from the corpus data
    doc_list = data["sequence"].tolist()
    return doc_list

#build_kmers returns the list of shingles of size 'ksize' from each document with id pasrsed to the parameter 'sequence'
def build_kmers(sequence, ksize):
    kmers = []
    n_kmers = len(sequence) - ksize + 1

    for i in range(n_kmers):
        kmer = sequence[i:i + ksize]
        kmers.append(kmer)

    return kmers

#set of unique shingles obtained from each document in the corpus
shingles={}

# shingle returns the set of shingles stored and dumped in a json file 'shingle_list.json'
def shingle():
    sequences = parsing_data()
    cnt=0
    shinglesInDocWords = set()
    PostingDict = {}
    for i in sequences:

        print("Shingling doc " + str(cnt+1))
        shingle_list=build_kmers(i,7)

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

    #print(shinglesInDocWords)
    list_of_unique_shingles=[]
    count=0
    for i in shinglesInDocWords:
        list_of_unique_shingles.insert(count,i)
        count+=1
    #print(count)
    # print(shinglesInDocWords)
    #with open("./shingle_list.json",'w') as ft:
    #    json.dump(list_of_unique_shingles, ft)
    #with open("./docwise_shingle_list.json",'w') as f:
    #    json.dump(docAsShingleSets, f)

    return docAsShingleSets, list_of_unique_shingles, PostingDict

'''
def invertedIndexMatrixGenerator(docsAsShingleSets, allShingles):
    
    #Paramters: docsAsShingleSets (The dictionary of documents with shingles)
    #           allShingles (all shingles generated till now from the corpus)
    #This function generates the posting list for each shingle in allShingles
    #It returns a dictionary of posting list for each shingle
    
    with open("./docwise_shingle_list.json") as data:
        docsAsShingleSets = json.load(data)
    with open("./shingle_list.json") as data:
        allShingles = json.load(data)
    print("Generating Inverted Index\n")
    invertedIndexTable = {}
    #allShingles = list(set(allShingles))
    cnt1 = 0
    for eachShingle in allShingles:
        postingsList = []
        cnt1 += 1
        print(cnt1)
        for j in docsAsShingleSets:
            if (eachShingle in docsAsShingleSets[j]): # If shingle in present in jth document,j is added to the list
                #try:
                postingsList.append(j)
                #except:
                    #postingsList = {j}
        invertedIndexTable[eachShingle] = postingsList # Inverted index table for each shingle is made
    print(invertedIndexTable)
    # with open("./inverted_table.json",'w') as ft:
    #     json.dump(invertedIndexTable, ft)
    
'''
def matrixGenerator(allShingles, invertedIndexTable):
    '''
    #Parameters: allShingles (list of all shingles in the corpus)
    #            invertedIndexTable (dictionary with posting list for each shingle)
    #This function indexes the shingles and returns a boolean matrix of shingle versus document
    '''
    postlist=[]
    index_matrix = {}
    index = 0
    # indexing the shingles
    print("Generating Boolean Matrix\n")
    for shingle in allShingles:
        index_matrix[shingle] = index
        index += 1

    # shingle document matrix
    matrix = np.zeros([len(allShingles), 1680], dtype=int)
    for shingle in allShingles:
        postlist = invertedIndexTable[shingle]
        for d in postlist:
            matrix[index_matrix[shingle]][int(d)] = 1  # Boolean value true for that document corresponding to a shingle

    return matrix

def pickRandomCoeffs(k, maxval):
    '''
    #Parameters: k (Number of random values)
    #            maxval (maximum value for randint)
    # This function returns k number of unique random values
    '''


  # Create a list of 'k' random values.
    randList = []

    while k > 0:
    # Get a random shingle ID.
        randIndex = random.randint(0, maxval)

    # Ensure that each random number is unique.
        while randIndex in randList:
            randIndex = random.randint(0, maxval)

    # Add the random number to the list.
        randList.append(randIndex)
        k = k - 1

    return randList

def find_sign_matrix(matrix, numOfShingles):
    '''
    #Parameters: matrix (boolean matrix of shingles vs docs)
    #            numOfShingles (total number of shingles in corpus)
    #This function picks two random coefficient values and genrates a hashfunction of form h(x)=(ax+b)%c
    #All values are initialised to infinities and each row is mapped to lowest hash function that has a
    #boolean true for that shingle. This new matrix called sigmatrix is returned.
    # example
    # matrix= [[1, 0, 0, 1],
    #         [0, 0, 1, 0],
    #         [0, 1, 0, 1],
    #         [1, 0, 1, 1],
    #         [0, 0, 1, 0]]
    # coeffA = [1,1]
    # coeffB = [1,3]
    # c = 5
    # required output is [[1, 3, 0, 1], [0, 2, 0, 0]]
    '''
    print("Generating signature Matrix\n")
    c = numOfShingles
    while not sympy.isprime(c):
        c += 1

    coeffA = pickRandomCoeffs(NUM_HASH_FUNCS, numOfShingles-1)
    coeffB = pickRandomCoeffs(NUM_HASH_FUNCS, numOfShingles-1)

    rows, cols, sigrows = len(matrix), len(matrix[0]), len(coeffA)
    # initialize signature matrix with maxint
    sigmatrix = np.full((sigrows, cols), INF)
    

    for r in range(rows):
        hashvalue = []
        for h in range(sigrows):
            hashvalue.append((coeffA[h] + coeffB[h]*r) % c)# Hash each row
        # if data != 0 and signature > hash value, replace signature with hash value
        for col in range(cols):
            if matrix[r][col] == 0:
                continue
            for i in range(sigrows):
                if sigmatrix[i, col] > hashvalue[i]:
                    sigmatrix[i, col] = hashvalue[i]
    print("Signature matrix\n")
    print(sigmatrix)
    return sigmatrix
#call the required functions

#parsing_data()
docsAsShingleSets, allShingles, PostingDict = shingle()
#print(PostingDict)
matrix = matrixGenerator(allShingles,PostingDict)
print(matrix)
sign_matrix = find_sign_matrix(matrix,len(allShingles))
#with open("./inverted_table.json",'w') as ft: 
#    json.dump(PostingDict, ft)
print("Dumped successfully")
#invertedIndexMatrixGenerator()

