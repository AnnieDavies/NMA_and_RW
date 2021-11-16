# Folder: PapakonstantinouAlgorithm

3 codes containing 3 versions of the evidence stream algorithm in Papakontaninou et al (2018); Shortest, Random and Average.

Codes: (A Davies 2021)

1. Algorithmm-Shortest.cpp

  This code was used to obtain the following results in the paper NMA and RW:
  
    a. 4th column of Table 3 (evidence streams)
    
    b. 3rd column of Table 4 (proportion contributions)
    
  
C++ code implementing the version of the algorithm that selects paths in order from shortest to longest.

2. Algorithmm-Random.cpp

  This code was used to obtain the following results in the paper NMA and RW:
  
    a. 5th column of Table 3 (evidence streams)
    
    b. 4th column of Table 4 (proportion contributions)
  
C++ code implementing the version of the algorithm that selects paths at random at each iteration. 

Each run of this code gives a different result. 

3. Algorithmm-Average.cpp

  This code was used to obtain the following results in the paper NMA and RW:
  
    a. 6th column of Table 3 (evidence streams)
    
    b. 5th column of Table 4 (proportion contributions)
  
  C++ code implementing the version of the algorithm that averages over the results of the Random algorithm for many realisations. 
  
  The integer MAX controls how many iterations of the algorithm are performed. 

Files:

1. Hflow13.txt

Describes the evidence flow through each edge from node 1 to node 3 on the evidence flow network for the Depression data set (Rucker and Schwarzer 2014)

Flow is obtained by the row of the hat matrix correpsonding to comparison 13 (see row 1 Hat_agg_netmeta.xlsx in RCodeDataResults)

(or, equivalently, using the analytical RW approach - see folder EvidenceFlow)

N = 11 (number of nodes)

The length of this file is N(N-1)/2.

Every possible pair of nodes is included.

If a certain pair of nodes are not connected by an edge in the evidence flow network then the flow is 0. 

