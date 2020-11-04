# DESIGN DOCUMENT

### **IMPLEMENTING LOCALITY SENSITIVE HASHING ON DNA DATASET WITH A QUERY AND FINDING THE SIMILAR DOCUMENTS**

### **SUBJECT- INFORMATION RETRIEVAL (CS F469)**

### **Group Members-**

RUPSA DHAR 2018A7PS0376H

KESHAV KABRA 2018AAPS0527H

MEGANAA REDDY 2017B3A70973H

### **Description of the System**

Our objective is to recommend the user with the similar sequences of DNA bases, from the corpus of DNA dataset with a threshold similarity given by the user. We exploit the concept of LSH , to find the documents which are closer to the query and display them to the user for various distance measures.

LSH refers to a family of functions to hash data points into buckets so that data points near each other are located in the same buckets with high probability, while data points far from each other are likely to be in different buckets. This makes it easier to identify observations with various degrees of similarity.

### **Data Structures Used**

shinglesInDocWords = set()#set to store unique hashed shingles of a doc,to avoid repetition

PostingDict = {} # Dictionary to store the hashed shingle mapped to a set of document id&#39;s which belong to it

docIDlist=set() #set of id&#39;s all the documents of a corpus

docAsShingleSets={} #docAsShingleSets is a dictionary of unique hashed\_shingles in a document with index i

list\_of\_unique\_shingles=[]#list of all unique shingles

invertedIndexTable (dictionary with posting list for each shingle)

matrix #boolean matrix with list\_of\_unique\_shingles as rows documents as columns, 1 indicates document contains shingle, otherwise 0

sigmatrix #signature matrix with no\_of\_rows=no\_of\_hash\_functions and no\_of\_columns=no\_of\_documents

### **Similarity Measure Used**
 **Jaccard coefficient measure** - It is a number between 0 and 1. It is defined as the number of elements in the intersection of two sets A and B divided by the number of elements in their union. The higher the coefficient, more is the similarity. (Lesser distance)
 
### **Runtime for Different Distance Measures**

Number of Documents = 1680, Number of Hash Functions = 100 , Shingles - 4 shingles

Threshold - 0.6

Optimized number of bands - 30

Jaccard - 0.2393sec

Total Runtime - 0.8393 sec


