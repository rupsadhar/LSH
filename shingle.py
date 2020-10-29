
import os
import json
import pandas as pd
import binascii as bin
data_list=[]


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
        shingle_list=build_kmers(i,4)

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
            if shingle_word not in PostingDict:
                PostingDict[shingle_word] = set()  

            PostingDict[shingle_word].add(cnt)

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

#call the required functions

#parsing_data()
docsAsShingleSets, allShingles, PostingDict = shingle()
print(PostingDict)
#with open("./inverted_table.json",'w') as ft: 
#    json.dump(PostingDict, ft)
print("Dumped successfully")
#invertedIndexMatrixGenerator()

