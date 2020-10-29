
import os
import json
import pandas as pd
data_list=[]
def parsing_data() :
    data = pd.read_csv("./dna_data/chimp_data-noN.txt", sep = "	")

    doc_list = data["sequence"].tolist()
    return doc_list

def build_kmers(sequence, ksize):
    kmers = []
    n_kmers = len(sequence) - ksize + 1

    for i in range(n_kmers):
        kmer = sequence[i:i + ksize]
        kmers.append(kmer)

    return kmers

shingles={}
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
    print(count)
    # print(shinglesInDocWords)
    with open("./shingle_list.json",'w') as ft:
        json.dump(list_of_unique_shingles, ft)


#parsing_data()
shingle()