# DESIGN DOCUMENT

### **IMPLEMENTING LOCALITY SENSITIVE HASHING ON DNA DATASET WITH A QUERY AND FINDING THE SIMILAR DOCUMENTS**

### **SUBJECT- INFORMATION RETRIEVAL (CS F469)**

### **Group Members-**

RUPSA DHAR 2018A7PS0376H

KESHAV KABRA 2018AAPS0527H

MEGANAA REDDY 2017B3A70973H

### **Description of the System**

Our objective is to recommend the user with the similar sequences of DNA bases, from the corpus of DNA dataset with a threshold similarity given by the user. We exploit the concept of LSH , to find the documents which are closer to the query and display them to the user for various distance measures.

### **Data Structures Used**

shinglesInDocWords = set()#set to store unique hashed shingles of a doc,to avoid repetition

PostingDict = {} # Dictionary to store the hashed shingle mapped to a set of document id&#39;s which belong to it

docIDlist=set() #set of id&#39;s all the documents of a corpus

docAsShingleSets={} #docAsShingleSets is a dictionary of unique hashed\_shingles in a document with index i

list\_of\_unique\_shingles=[]#list of all unique shingles

invertedIndexTable (dictionary with posting list for each shingle)

matrix #boolean matrix with list\_of\_unique\_shingles as rows documents as columns, 1 indicates document contains shingle, otherwise 0

sigmatrix #signature matrix with no\_of\_rows=no\_of\_hash\_functions and no\_of\_columns=no\_of\_documents

### **Functions Used**
1.**parsing_data(inputQuery)**
Parses data form txt file and converts it to a list.
:return: doc_list: list of documents list from the dataset

### Note
* Number of Documents = 1680
* Number of Hash Functions = 100 
* Shingle size = 7
* Threshold = 0.6
* Optimized number of bands = 18

### **Runtime for Different Distance Measures**
* Shingling =  60.93597 s
* Generating Signature matrix = 113.13688 s
* LSH = 0.18276 s
* Jaccard Similarity =  0.000587 s


### **Similarity Measure Used**
 **Jaccard coefficient measure** - It is a number between 0 and 1. It is defined as the number of elements in the intersection of two sets A and B divided by the number of elements in their union. The higher the coefficient, more is the similarity. (Lesser distance)
