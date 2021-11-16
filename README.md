# Network meta-analysis and random walks: Data, codes and results
Data, codes and results for the manuscript Network meta-analysis and random walks.

The files are separated into five folders each with their own more detailed README file. 
Each code is used to obtain a set of results presented in the manuscript NMA and RW. 
The codes are applied to the Depression data set described in: 
Rücker G, Schwarzer G. Reduce dimension or reduce weights? Comparing two approaches to multi-arm studies in networkmeta-analysis.Stat Med2014; 33(25): 4353-4369

The relevant data is given in files in the appropriate folders. 

All codes were written by A Davies (2021).

## Folders:

## 1. EvidenceFlow

Contains two codes for working out evidence flow networks. 

  a. AggregateHatMatrix.cpp
  
       Performs a frequentist NMA using the graph theoretical approach both as a one-step and two-step method (method described in Section 3). 
       
       Yields the hat matrix of the aggregate model - each row describes an evidence flow network for the comparison it represents.
       
       This code was used to obtain the following results in the paper NMA and RW:
       
         i. The aggregate weight matrix in Appendix F (= aggregate weights listed in Figure 6)
         
         ii. The hat matrix of the aggregate model listed in Appendix F
         
 b. EvidenceFlow_RWanalytical.cpp
 
       Works out the flow of evidence for the comparison source-sink using the analytical RW approach (method described in Section 4.3 and Appendix E).
       
       This code was used to obtain the following results in the paper NMA and RW:
       
         i. The transition matrix for a RW on the aggregate network (Eq (26))
         
         ii. The evidence flow network for comparison 1-3 of the Depression data shown in Figure 7
         
 c. Files
 
      N = 11 (number of nodes)
      
      T = 26 (number of trials)
      
      K = 20 (number of edges)
      
      Q = 47 (total number of pairwise comparisons)
      
      Data files provided are for the Depression data described in Rucker and Schwarzer (2014). 
      
        i. M.txt = list (length T) of the number of treatments in each trial
        
        ii. x.txt = list (length Q) of the relative treatment effects (as log odds ratios) measured in each trial
        
        iii. SD.txt = list (length Q) of the standard deviation associated with each value of x
        
        iv. B.txt = the edge-incidence matrix (QxN) of the network of treatments and trials, describes which treatments are in each trial
        
        v. Ba.txt = the edge-incidence matrix (KxN) of the aggregate network, describes which treatment comparisons are represented by each edge
        

## 2. EvidenceStreams

  Obtaining the evidence streams and proportion contributions from the RW approach.
  
  a. Streams_contributions_RWanalytical.cpp 
  
      Code to obtain the evidence streams and proportion contributions using the analytical RW approach (method described in Section 5.3).
      
      This code was used to obtain the following results in the paper NMA and RW:
      
        i. The transition matrix on the evidence flow network (Eq 27)
        
        ii. Columns 2 and 3 of Table 3 (evidence streams)
        
        iii. Column 2 of Table 4 (proportion contribution)
        
  b. Hflow13.txt
  
      Describes the evidence flow through each edge from node 1 to node 3 on the evidence flow network for the Depression data set.
      
      Flow is obtained by the row of the hat matrix correpsonding to comparison 13
      
      (or, equivalently, using the analytical RW approach - see folder EvidenceFlow)
      
      The length of this file is N(N-1)/2. Every possible pair of nodes is included.
      
      If a certain pair of nodes are not connected by an edge in the evidence flow network then the flow is 0. 
      
      
## 3. PapakonstantinouAlgorithm

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

## 4. RCodeDataResults

R code, data and results to obtain some results in the paper using the R package netmeta. 

  a. HatAgg_Contrib_Depression.R

       R code that reads in the Depression data and works out the following using netmeta:
            i.   Aggregate hat matrix (hat matrix in Section H.1 of the Supplementary Material)
            ii.  Contribution matrix using the random walk approach (1st row gives column 2 in Table 4)
            iii. Contribution matrix using the `Shortest' algorithm approach (1st row gives column 3 in Table 4)

  b. Rucker2014.txt

       Data Frame for the Depression data set (used in Rucker and Schwarzer (2014)).  

  c. Results from the code HatAgg_Contrib_Depression.R saved as excel files.
            i.   Hat_agg_netmeta.xlsx (aggregate hat matrix for a fixed effect model)
            ii.  Contrib_randomwalk_netmeta.xlsx (contrbution matrix for a fixed effect model using the random walk approach)
            iii. Contrib_shortestpath_netmeta.xlsx (contrbution matrix for a fixed effect model using the `Shortest' algorithm)

## 5. Results

Files containing the results from the paper. 

Some of these are presented in the paper but here they are quoted to higher precision. 

Others are not explicitly presented in the paper.

  a. Algorithm_Results
  
          Results obtained from the Algorithm in Papakonstantinou et al (2018).
          
          Files refer to the following results in the paper NMA and RW
          
            i. Av-Prop.txt = 5th column in Table 4 = proportion contributions as estimated by the algorithm 'Average' (code: Algorithmm-Average.cpp)
            
            ii. Av-Stream.txt = 6th column in Table 3 = streams as estimated by the algorithm 'Average'  (code: Algorithmm-Average.cpp)
          
            iii. Short-Prop.txt = 3rd column in Table 4 = proportion contributions as estimated by the algorithm 'Shortest' (code: Algorithmm-Shortest.cpp)
            
            iv. Short-Stream.txt = 4th column in Table 3 = streams as estimated by the algorithm 'Shortest'  (code: Algorithmm-Shortest.cpp), paths are labelled by their composite nodes
            
  b. NMA_results
  
          Results (and intermediate steps) from the two-step graph theoretical NMA model (as described in Section 3).
          
          All produced from the code AggregateHMatrix.cpp
          
            i. Hat.txt = hat matrix of the aggregate model (calculated from Eq (5), presented in Appendix F)
            
            ii. LOR_dir.txt = theta^(dir) = direct estimates as log odds ratios (calculated from Eq (1))
            
            iii. LOR_net.txt = theta^(net) = network estimates as log odds ratios (calculated from Eq (4))
            
            iv. OR_dir.txt = theta^(dir) = direct estimates as odds ratios (exponential of LOR_dir.txt)
            
            v. OR_net.txt = theta^(net) = network estimates as odds ratios (exponential of LOR_net.txt)
            
            vi. Wagg.txt = aggregate weight matrix (calculated from Eq (2), presented in Appendix F)
            
            vii. Wfull.txt = full adjusted weight matrix (calculated from the adjustment method for multi-arm trials - see Rucker and Schwarzer 2014)
  
  c. RW_results

          Results from the random walk approach.
          
            i. RW_flow_13.txt = the evidence flow for comparison 1-3 of the Depression data as estimated by the analytical RW approach (see Section 4.3 and Appendix E, presented in Figure 7, code: EvidenceFlow_RWanalytical.cpp)
            
            ii. RW_prop.txt = proportion contributions to the network comparison 1-3 in the Depression data as estimated by the RW approach (see Section 5.3, presented in 2nd column of Table 4, code: Streams_contributions_RWanalytical.cpp)
            
            iii. RW_stream.txt = evidence streams (flow through paths) for the network comparison 1-3 in the Depression data set as estimated by the RW approach (see section 5.3, presented in 3rd column of Table 3, code: Streams_contributions_RWanalytical.cpp)
            
            iv. T_13.txt = The transition matrix for a RW on the aggregate network for the Depression data (Figure 6) moving from node 1 to node 3 (calculated from Eq (12) and Wagg.txt, presented in Eq (26), code: EvidenceFlow_RWanalytical.cpp)
            
            v. U_13.txt = The transition matrix for a RW on the evidence flow network for comaprison 1-3 in the Depression data (Figure 7) moving from node 1 to node 3 (calculated from Eq (20) and RW_flow_13.txt, presented in Eq (27), code: Streams_contributions_RWanalytical.cpp)
