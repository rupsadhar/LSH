
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
    for i in sequences:
        print("Shingling doc " + str(cnt+1))
        shingles[cnt]=build_kmers(i,10)
        cnt+=1
    with open("./shingles.json",'w') as ft:
        json.dump(shingles, ft)


shingle()