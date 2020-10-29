
import os
import json
import pandas as pd
data_list=[]
def parsing_data() :
    # f=open("human_data.txt","r")
    # lines=f.readlines()
    # result=[]
    # for x in lines:
    #     result.append(x.split('\t')[0])
    # f.close()
    # with open("doc_as_strings.json",'w') as f1:
    #     json.dump(result,f1)
    
    data = pd.read_csv("./dna_data/human_data.txt", sep = "	")

    doc_list = data["sequence"].tolist()
# print(doc_list)

    with open("./dna_data/human_data.json",'w') as ft:
        json.dump(doc_list, ft)

def build_kmers(sequence, ksize):
    kmers = []
    n_kmers = len(sequence) - ksize + 1

    for i in range(n_kmers):
        kmer = sequence[i:i + ksize]
        kmers.append(kmer)

    return kmers

shingles={}
def shingle():
    with open("./dna_data/human_data.json") as data:
        sequences = json.load(data)
    cnt=0
    shinglesInDocWords = set()
    for i in sequences:
        shingle_list=build_kmers(i,1)
        for shingle_word in shingle_list:
            if shingle_word not in shinglesInDocWords:
                shinglesInDocWords.add(shingle_word)
                cnt+=1
    print(cnt)
    list_of_unique_shingles=[]
    count=0
    for i in shinglesInDocWords:
        list_of_unique_shingles.insert(count,i)
        count+=1
    # print(shinglesInDocWords)
    with open("./shingle_list.json",'w') as ft:
        json.dump(list_of_unique_shingles, ft)

parsing_data()
shingle()