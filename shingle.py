
import os
import json
import pandas as pd
import binascii as bin
import numpy as np
import sympy as sympy
import random as random
import math
data_list=[]
INF = 2**32
NUM_HASH_FUNCS=100
threshold = 0.6
B_ROWS = 5
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
docIDlist=set()
# shingle returns the set of shingles stored and dumped in a json file 'shingle_list.json'
def shingle():
    sequences = parsing_data()
    cnt=0
    shinglesInDocWords = set()
    PostingDict = {}
    for i in sequences:

        print("Shingling doc " + str(cnt+1))
        docIDlist.add(cnt+1)
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

    return docAsShingleSets, list_of_unique_shingles, PostingDict, docIDlist

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
    check=0
    matrix = np.zeros([len(allShingles), 1680], dtype=int)
    for shingle in allShingles:
        postlist = invertedIndexTable[shingle]
        for d in postlist:
            matrix[index_matrix[shingle]][int(d)] = 1  # Boolean value true for that document corresponding to a shingle
            check+=1
    print(check)
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
def getbestb(threshold,NUM_HASH_FUNCS, eps=1e0):
    '''
    #Parameters: threshold (difined threshold)
    #            NUM_HASH_FUNCS (number of hash functions)
    #            eps
    # Returns the best value for b by solving an equation
    '''
    for b in range(1, NUM_HASH_FUNCS+1):
        opt = b*math.log10(b)
        val = -1 * NUM_HASH_FUNCS * math.log10(threshold)
        if opt > val-eps and opt < val+eps:
            print("Using number of bands : %d" % (np.round(b)))
            return np.round(b)
def lsh(B_BANDS, docIdList, sig):
    '''
    #Parameters: B_BANDS (Number of bands in signature matrix)
    #            docIdList (List of document ids)
    #            sig (signature matrix)
    #This function first divides the signature matrix into bands and hashes each column onto buckets.
    #This hashing is called Locality Sensitive Hashing.
    #This function returns the list of document to its hash along with the buckets
    '''
    n = NUM_HASH_FUNCS
    b = getbestb(threshold,NUM_HASH_FUNCS)
    r = n / b

    d = 1680
    # Array of dictionaries, each dictionary is for each band which will hold buckets for hashed vectors in that band
    buckets = np.full(b, {})
    # Mapping from docid to h to find the buckets in which document with docid was hashed
    docth = np.zeros((d, b), dtype=int)  # doc to hash
    for i in range(b):
        for j in range(d):
            low = int(i*r) # First row in a band
            high = min(int((i+1)*r), n)# Last row in current band
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
    '''
    Jackard similarity
    '''
    doc1 = sign_matrix[:,doc1]
    doc2 = sign_matrix[:,doc2]
    intersection = sum(bool(x) for x in np.logical_and(doc1, doc2))
    return (intersection / len(doc1))


def euclidean_distance(x, y, r=2.0):
    '''
    Euclidean distance
    '''
    A = np.linalg.norm(x)
    B = np.linalg.norm(y)
    x = np.divide(x, A)
    y = np.divide(y, B)
    try:
         return np.power(np.sum(np.power(np.subtract(x, y), r)), 1.0/r)
    except (ValueError,ZeroDivisionError):
         print('Please, enter only even values for "r > 0".')
    except IndexError:
         print('Please, the sets must have the same size.')

def cosine_distance(x,y):
    '''
    Cosine Distance
    '''
    prodAB = np.dot(x,y)
    #zeros = np.zeros(len(x))
    A = np.linalg.norm(x)
    B = np.linalg.norm(y)
    return prodAB / (A*B)


def get_similarcos(dn,docIdList,buckets,docth,docsAsShingleSets,sign_matrix):
    '''
    Similarity for cosine distance
    '''
    if dn not in docIdList:
        raise KeyError('No document with the given name found in the corpus.')

    docid = dn
    # Collection of documents similar to docid
    c = []
    # taking union of all buckets in which docid is present
    for b, h in enumerate(docth[docid]):
        c.extend(buckets[b][h])
    c = set(c)
    print(c)

    # Similar documents
    sim_list = []
    for doc in c:
        if doc == docid:
            continue
        sim = cosine_distance(sign_matrix[:,dn],
                              sign_matrix[:,doc])
        sim_list.append((sim, doc))
    sim_list.sort()
    return sim_list

def get_similar(dn,docIdList,buckets,docth,docsAsShingleSets,sign_matrix):
    '''
    #Parameters: dn (The query document number)
    #            docIdList (List of doc ids)
    #            buckets (List of buckets)
    #            docth (doc to hash list)
    #            docAsShingleSets
    #            Signature Matrix   
    # This function finds similar documents given a query document after hashing and bucketing the query document
    # It also evaluates based on various similarity criterion, namely, Jacard similarity, Euclidean distance
    # and cosine similarity
    # It returns a list of similar documents based on decreasing similarity amount
    '''
    if dn not in docIdList:
        raise KeyError('No document with the given name found in the corpus.')

    docid = dn
    # Collection of documents similar to docid
    c = []
    # taking union of all buckets in which docid is present
    for b, h in enumerate(docth[docid]):
        c.extend(buckets[b][h])
    c = set(c)
    print(c)

    # Similar documents
    sim_list = []
    for doc in c:
        if doc == docid:
            continue
        sim = jacsim(docid, doc, docsAsShingleSets,sign_matrix)
        sim_list.append((sim, doc))
    sim_list.sort(reverse=True)
    return sim_list

def get_similareucdis(dn,docIdList,buckets,docth,docsAsShingleSets,sign_matrix):
    '''Similarity For Euclidean Distance'''
    if dn not in docIdList:
        raise KeyError('No document with the given name found in the corpus.')

    docid = dn
    # Collection of documents similar to docid
    c = []
    # taking union of all buckets in which docid is present
    for b, h in enumerate(docth[docid]):
        c.extend(buckets[b][h])
    c = set(c)
    print(c)

    # Similar documents
    sim_list = []
    for doc in c:
        if doc == docid:
            continue
        sim = euclidean_distance(sign_matrix[:,dn],sign_matrix[:,doc])
        sim_list.append((sim, doc))
    sim_list.sort()
    return sim_list
#call the required functions

#parsing_data()
docsAsShingleSets, allShingles, PostingDict, docIDlist = shingle()
#print(PostingDict)
matrix = matrixGenerator(allShingles,PostingDict)
print(matrix)
sign_matrix = find_sign_matrix(matrix,len(allShingles))
BANDS=20
docth,buckets = lsh(BANDS,docIDlist,sign_matrix)

print(docIDlist)
inputDocID = input("enter the doc ID you want to know similarities of : ")



#Using Jaccard Similarity
sim_docs = get_similar(int(inputDocID),docIDlist,buckets,docth,docsAsShingleSets,sign_matrix)

print("Calculating Jaccard similarities....\n")

found = 0
for sim, doc in sim_docs:
    if sim >= threshold:
        found = 1
        print('Document Name: ' + str(doc), 'Similarity: ' + str(sim) + '\n')
if found == 0:
    print("NO similar docs for the given threshold")
    
sim_docs1 = get_similarcos(int(inputDocID),docIDlist,buckets,docth,docsAsShingleSets,sign_matrix)

print("Calculating Cosine similarities....\n")

found = 0
for sim, doc in sim_docs1:
    if sim >= threshold:
        found = 1
        print('Document Name: ' + str(doc), 'Similarity: ' + str(sim) + '\n')
if found == 0:
    print("NO similar docs for the given threshold")

sim_docs2 = get_similareucdis(int(inputDocID),docIDlist,buckets,docth,docsAsShingleSets,sign_matrix)
    
print("Calculating Euclidean Distance....\n")

found = 0
t = np.sqrt(2 - 2*threshold)
for sim, doc in sim_docs2:
    if sim < t:
        found = 1
        print('Document Name: ' + str(doc), 'Similarity: ' + str(sim) + '\n')
if found == 0:
    print("NO similar docs for the given threshold")   


       
#with open("./inverted_table.json",'w') as ft: 
#    json.dump(PostingDict, ft)
print("Dumped successfully")
#invertedIndexMatrixGenerator()

