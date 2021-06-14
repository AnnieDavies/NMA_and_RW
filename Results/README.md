# Results

Files containing the results from the paper. 

Some of these are presented in the paper but here we quote them to a higher precision. 

Others are not explicitly presented in the paper.

1. Algorithm_Results
  
          Results obtained from the Algorithm in Papakonstantinou et al (2018).
          
          Files refer to the following results in the paper NMA and RW
          
            i. Av-Prop.txt = 5th column in Table 4 = proportion contributions as estimated by the algorithm 'Average' (code: Algorithmm-Average.cpp)
            
            ii. Av-Stream.txt = 6th column in Table 3 = streams as estimated by the algorithm 'Average'  (code: Algorithmm-Average.cpp)
          
            iii. Short-Prop.txt = 3rd column in Table 4 = proportion contributions as estimated by the algorithm 'Shortest' (code: Algorithmm-Shortest.cpp)
            
            iv. Short-Stream.txt = 4th column in Table 3 = streams as estimated by the algorithm 'Shortest'  (code: Algorithmm-Shortest.cpp), paths are labelled by their composite nodes
            
2. NMA_results
  
          Results (and intermediate steps) from the two-step graph theoretical NMA model (as described in Section 3.1).
          
          All produced from the code AggregateHMatrix.cpp
          
            i. Hat.txt = hat matrix of the aggregate model (calculated from Eq (5), presented in Appendix F)
            
            ii. LOR_dir.txt = theta^(dir) = direct estimates as log odds ratios (calculated from Eq (1))
            
            iii. LOR_net.txt = theta^(net) = network estimates as log odds ratios (calculated from Eq (4))
            
            iv. OR_dir.txt = theta^(dir) = direct estimates as odds ratios (exponential of LOR_dir.txt)
            
            v. OR_net.txt = theta^(net) = network estimates as odds ratios (exponential of LOR_net.txt)
            
            vi. Wagg.txt = aggregate weight matrix (calculated from Eq (2), presented in Appendix F)
            
            vii. Wfull.txt = full adjusted weight matrix (calculated from the adjustment method for multi-arm trials - see Rucker and Schwarzer 2014)
 
  
3. RW_results

          Results from the random walk approach.
          
            i. RW_flow_13.txt = the evidence flow for comparison 1-3 of the Depression data as estimated by the analytical RW approach (see Section 4.3 and Appendix E, presented in Figure 7, code: EvidenceFlow_RWanalytical.cpp)
            
            ii. RW_prop.txt = proportion contributions to the network comparison 1-3 in the Depression data as estimated by the RW approach (see Section 5.3, presented in 2nd column of Table 4, code: Streams_contributions_RWanalytical.cpp)
            
            iii. RW_stream.txt = evidence streams (flow through paths) for the network comparison 1-3 in the Depression data set as estimated by the RW approach (see section 5.3, presented in 3rd column of Table 3, code: Streams_contributions_RWanalytical.cpp)
            
            iv. T_13.txt = The transition matrix for a RW on the aggregate network for the Depression data (Figure 6) moving from node 1 to node 3 (calculated from Eq (12) and Wagg.txt, presented in Eq (26), code: EvidenceFlow_RWanalytical.cpp)
            
            v. U_13.txt = The transition matrix for a RW on the evidence flow network for comaprison 1-3 in the Depression data (Figure 7) moving from node 1 to node 3 (calculated from Eq (20) and RW_flow_13.txt, presented in Eq (27), code: Streams_contributions_RWanalytical.cpp)


Order of pairwise comparisons in Hat.txt (rows and columns), LOR_dir.txt, LOR_net.txt, OR_dir.txt. OR_net.txt, Wagg.txt (rows and columns), RW_flow_13.txt:
  
  (1→3,		1→6,		1→7,		1→9,		1→11,		2→6,		2→8,		2→11,		3→4,		3→5,		3→6, 		3→9,		4→9,		5→9,		6→7,		6→8,		6→9,		6→11,		7→9,	7→10)
  
Order of pairwise comparisons in Av-Prop.txt, Short-Prop.txt, RW_prop.txt:

  (1→2, 1→3, 1→4, 1→5, 1→6, 1→7, 1→8, 1→9, 1→10, 1→11, 2→3, 2→4, 2→5, 2→6, 2→7, 2→8, 2→9, 2→10, 2→11, 3→4, 3→5, 3→6, 3→7, 3→8, 3→9, 3→10, 3→11, 4→5, 4→6, 4→7, 4→8, 4→9, 4→10, 4→11, 5→6, 5→7, 5→8, 5→9, 5→10, 5→11, 6→7, 6→8, 6→9, 6→10, 6→11, 7→8, 7→9, 7→10, 7→11, 8→9, 8→10, 8→11, 9→10, 9→11, 10→11)

Order of paths in Av-Stream.txt and RW-Stream.txt as in Table 3 (NB: in Short-Stream.txt order is different and paths are labelled accordingly in the file):

1→3,   
1→6→3,   
1→6→9→3,   
1→6→9→4→3,   
1→6→9→5→3,   
1→7→6→3,   
1→7→6→9→3,   
1→7→6→9→4→3,   
1→7→6→9→5→3,   
1→7→9→3,   
1→7→9→4→3,   
1→7→9→5→3,   
1→9→3,   
1→9→4→3,   
1→9→5→3,   
1→11→2→6→3,   
1→11→2→6→9→3,   
1→11→2→6→9→4→3,   
1→11→2→6→9→5→3,   
1→11→2→8→6→3,  
1→11→2→8→6→9→3,   
1→11→2→8→6→9→4→3,   
1→11→2→8→6→9→5→3,   
1→11→6→3,   
1→11→6→9→3,   
1→11→6→9→4→3,  
1→11→6→9→5→3,   




  
