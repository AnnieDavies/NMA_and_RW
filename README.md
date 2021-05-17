# NMA_and_RW
Codes for the manuscript Network meta-analysis and random walks.

The codes are separated into three files each with their own more detailed README file. 
Each code is used to obtain a set of results presented in the manuscript NMA and RW. 
The codes are applied to the Depression data set described in: 
Rücker G, Schwarzer G. Reduce dimension or reduce weights? Comparing two approaches to multi-arm studies in networkmeta-analysis.Stat Med2014; 33(25): 4353-4369

The relevant data is given in files in the appropriate folders. 

All codes are in C++ and were written by A Davies (2021).

Folders:

1. EvidenceFlow

Contains two codes for working out evidence flow networks. 

  a. AggregateHatMatrix.cpp
       Performs a frequentist NMA using the graph theoretical approach both as a one-step and two-step method (method described in Section 3.1). 
       
       Yields the hat matrix of the aggregate model - each row describes an evidence flow network for the comparison it represents.
       
       This code was used to obtain the following results in the paper NMA and RW:
       
         i. The aggregate weight matrix in Appendix D (= aggregate weights listed in Figure 6)
         
         ii. The hat matrix of the aggregate model listed in Appendix D
         
 b. EvidenceFlow_RWanalytical.cpp
 
       Works out the flow of evidence for the comparison source-sink using the analytical RW approach (method described in Section 3.4 and Appendix C.2).
       
       This code was used to obtain the following results in the paper NMA and RW:
       
         i. The transition matrix for a RW on the aggregate network (Eq (23))
         
         ii. The evidence flow network for comparison 1-3 of the Depression data shown in Figure 7
         
 c. Files
 
      N = 11 (number of nodes)
      
      T = 26 (number of trials)
      
      K = 20 (number of edges)
      
      Q = 47 (total number of pairwise comparisons)
      
      Data files provided are for the Depression data described in Rucker and Schwarzer (2014). 
      
        i. M.txt = list (length T) of the number of treatments in each trial
        
        ii. x.txt = list (length Q) of the relative treatment effects measured in each trial
        
        iii. SD.txt = list (length Q) of the standard deviation associated with each value of x
        
        iv. B.txt = the edge-incidence matrix (QxN) of the network of treatments and trials, describes which treatments are in each trial
        
        v. Ba.txt = the edge-incidence matrix (KxN) of the aggregate network, describes which treatment comparisons are represented by each edge
        

2. EvidenceStreams

  Obtaining the evidence streams and proportion contributions from the RW approach.
  
  a. Streams_contributions_RWanalytical.cpp 
  
      Code to obtain the evidence streams and proportion contributions using the analytical RW approach (method described in Section 4.2).
      
      This code was used to obtain the following results in the paper NMA and RW:
      
        i. The transition matrix on the evidence flow network (Eq 24)
        
        ii. Columns 2 and 3 of Table 3 (evidence streams)
        
        iii. Column 2 of Table 4 (proportion contribution)
        
  b. Hflow13.txt
  
      Describes the evidence flow through each edge from node 1 to node 3 on the evidence flow network for the Depression data set.
      
      Flow is obtained by the row of the hat matrix correpsonding to comparison 13
      
      (or, equivalently, using the analytical RW approach - see folder EvidenceFlow)
      
      The length of this file is N(N-1)/2. Every possible pair of nodes is included.
      
      If a certain pair of nodes are not connected by an edge in the evidence flow network then the flow is 0. 
      
      
3. PapakonstaninouAlgorithm

  Obtaining the evidence streams and proportion contributions from the algorithm in:
  
  Papakonstantinou T, Nikolakopoulou A, Rücker G, et al. Estimating the contribution of studies in network meta-analysis:paths, flows and streams..F1000Research2018; 7: 610
  
  a. Algorithmm-Shortest.cpp
  
      The version of the algorithm that selects paths in order from shortest to longest.
      
      This code was used to obtain the following results in the paper NMA and RW:
      
        i. 4th column of Table 3 (evidence streams)
        
        ii. 3rd column of Table 4 (proportion contributions)
        
  b. Algorithm-Random.cpp
  
      The version of the algorithm that selects paths at random at each iteration.
      
      This code was used to obtain the following results in the paper NMA and RW:
      
        i. 5th column of Table 3 (evidence streams)
        
        ii. 4th column of Table 4 (proportion contributions)
        
  c. Algorithmm-Average.cpp
  
      The version of the algorithm that averages over the results of the Random algorithm for many realisations.
      
      This code was used to obtain the following results in the paper NMA and RW:
      
        i. 6th column of Table 3 (evidence streams)
        
        ii. 5th column of Table 4 (proportion contributions)  
        
  d. Hflow13.txt
  
      Describes the evidence flow through each edge from node 1 to node 3 on the evidence flow network for the Depression data set.
      
      Flow is obtained by the row of the hat matrix correpsonding to comparison 13
      
      (or, equivalently, using the analytical RW approach - see folder EvidenceFlow)
      
      The length of this file is N(N-1)/2. Every possible pair of nodes is included.
      
      If a certain pair of nodes are not connected by an edge in the evidence flow network then the flow is 0.
