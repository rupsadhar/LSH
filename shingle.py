
import os
import json
import pandas as pd
data_list=[]

#parsing_data() returns the list of documents list from the dataset
def parsing_data() :
    #read data as .csv 
    data = pd.read_csv("./dna_data/human_data.txt", sep = "	")
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
    for i in sequences:
        print("Shingling doc " + str(cnt+1))
        shingle_list=build_kmers(i,1)
        for shingle_word in shingle_list:
            if shingle_word not in shinglesInDocWords:
                shinglesInDocWords.add(shingle_word)
        cnt+=1
    print(cnt)

    print(shinglesInDocWords)
    list_of_unique_shingles=[]
    count=0
    for i in shinglesInDocWords:
        list_of_unique_shingles.insert(count,i)
        count+=1
    # print(shinglesInDocWords)
    with open("./shingle_list.json",'w') as ft:
        json.dump(list_of_unique_shingles, ft)

#call the required functions
parsing_data()
shingle()
