import pandas as pd 
import json
data = pd.read_csv("./dna_data/human_data.txt", sep = "	")

doc_list = data["sequence"].tolist()
# print(doc_list)

with open("./dna_data/human_data.json",'w') as ft:
     json.dump(doc_list, ft)