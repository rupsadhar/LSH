DESIGN DOCUMENT

IMPLEMENTING LOCALITY SENSITIVE HASHING ON DNA DATASET WITH A QUERY AND FINDING THE SIMILAR DOCUMENTS
SUBJECT- INFORMATION RETRIEVAL (CS F469)

Group Members-

•	RUPSA DHAR		2018A7PS0376H
•	KESHAV KABRA		2018AAPS0527H
•	MEGANAA REDDY 	2017B3A70973H

Description of the System

Our objective is to recommend the user with the similar sequences of DNA bases, from the corpus of DNA dataset with a threshold similarity given by the user. We exploit the concept of LSH , to find the documents which are closer to the query and display them to the user for various distance measures.

The Program can be broken into various subparts

	Various function definitions from the file shingle.py does the following 
•	parsing_data() : Parsing the dataset into set of documents separated by a line from the file ‘chimp_data_noN.txt’ from the folder dna_data inside the main directory
•	build_kmers(): returns the list of shingles of size 'ksize' from each document with id parsed to the parameter 'sequence'
•	shingle(): returns list of unique shingles, set of shingles within a document of index i, dictionary to store hashed_shingles in a document and a list of documents
•	matrixGenerator(): This function indexes the shingles and returns a boolean matrix of shingle versus document
•	pickRandomCoeffs():This function returns k number of unique random values
•	find_sign_matrix(): This function picks two random coefficient values and generates a hashfunction of form h(x)=(ax+b)%c
•	lsh(): This function returns the list of document to its hash along with the buckets

File to be Executed
	
$ python 3 shingle.py
Data Structures Used 
•	shinglesInDocWords = set() -> set to store unique hashed shingles of a doc,to avoid repetition

•	PostingDict = { } ->Dictionary to store the hashed shingle mapped to a set of document id's which belong to it
•	docIDlist=set() -> set of id's all the documents of a corpus
•	docAsShingleSets={ } -> docAsShingleSets is a dictionary of unique hashed_shingles in a document with index i
•	list_of_unique_shingles=[ ] -> ist of all unique shingles
•	invertedIndexTable (dictionary with posting list for each shingle)
•	matrix  -> boolean matrix with list_of_unique_shingles as rows documents as columns, 1 indicates document contains shingle, otherwise 0
•	sigmatrix -> signature matrix with no_of_rows=no_of_hash_functions and no_of_columns=no_of_documents
   
Machine specs:
1.	Processor: i5-8250U
2.	Ram: 8 GB DDR3
3.	OS: Windows 10 WSL2 kernel

Runtime for Different Distance Measures

Number of Documents = 1680, Number of Hash Functions = 100 , Shingles - 4 shingles 
Threshold - 0.6
Optimized number of bands - 30
Jaccard - 0.2393sec
Cosine -0.2094sec
Euclidean -0.0069sec
Total Runtime - 0.8sec

Distances Used

1.	Euclidean distance- A simple vector space euclidean distance is used here.
2.	Cosine distance-  Useful to remove the effects of size of document in determining the similarity level. The lower the angle between query and doc or higher the cosine coefficient between the vectors implies ,higher is the similarity between them. 
3.	Jaccard coefficient measure- It is a number between 0 and 1. It is defined as the number of elements in the intersection of two sets A and B divided by the number of elements in their union. The higher the coefficient, more is the similarity. (Lesser distance)





