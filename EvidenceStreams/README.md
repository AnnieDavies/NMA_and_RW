# Folder: EvidenceStreams
Code: (A Davies 2021)

1. Streams_contributions_RWanalytical.cpp

  Code to obtain the evidence streams and proportion contributions using the RW approach.
  
  This code was used to obtain the following results in the paper NMA and RW:
  
    a. The transition matrix on the evidence flow network (Eq 27)
    
    b. Columns 2 and 3 of Table 3 (evidence streams)
    
    c. Column 2 of Table 4 (proportion contribution)

File:

1. Hflow13.txt

  Describes the evidence flow through each edge from node 1 to node 3 on the evidence flow network for the Depression data set (Rucker and Schwarzer 2014)
  
  Flow is obtained by the row of the hat matrix correpsonding to comparison 13 (see row 1 of Hat_agg_netmeta.xlsx in RCodeDataResults). 
  
  (or, equivalently, using the analytical RW approach - see folder EvidenceFlow)
  
  N = 11 (number of nodes)
  
  The length of this file is N(N-1)/2.
  
  Every possible pair of nodes is included.
  
  If a certain pair of nodes are not connected by an edge in the evidence flow network then the flow is 0. 
