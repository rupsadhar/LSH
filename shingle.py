
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
    
    data = pd.read_csv("./human_data.txt", sep = "	")

    doc_list = data["sequence"].tolist()
# print(doc_list)

    with open("./human_data.json",'w') as ft:
        json.dump(doc_list, ft)

parsing_data()