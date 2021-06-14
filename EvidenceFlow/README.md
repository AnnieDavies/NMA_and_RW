# Folder: EvidenceFlow

N = number of nodes

T = number of trials

K = number of edges

Q = total number of pairwise comparisons



Codes: (written by A Davies 2021)
1. AggregateHatMatrix.cpp

  This code was used to obtain the following results in the paper NMA and RW:

    a. The aggregate weight matrix in Appendix F (= aggregate weights listed in Figure 6)

    b. The hat matrix of the aggregate model listed in Appendix F

  C++ code that performs a frequentist NMA using the graph theoretical approach
  
  Both the one-step and two-step versions are provided
  
  The following data files must be provided by the user
  
  M.txt = list (length T) of the number of treatments in each trial
  
  x.txt = list (length Q) of the relative treatment effects measured in each trial
  
  SD.txt = list (length Q) of the standard deviation associated with each value of x
  
  B.txt = the edge-incidence matrix (QxN) of the network of treatments and trials, describes which treatments are in each trial
  
  Ba.txt = the edge-incidence matrix (KxN) of the aggregate network, describes which treatment comparisons are represented by each edge


  The hat matrix of the aggregate model is calculate (Ha). 
  
  Each row of this matrix describes an evidence flow network for the comparison it represents. 


  The code can also be used to obtain the network estimates of the relative treatment effects. 
  
  From the one-step model these are 'expTheta' and from the two-step (aggregate) model these are 'expTheta_a'.
  
  The models are equivalent and therefore expTheta = expTheta_a
  

2. EvidenceFlow_RWanalytical.cpp

  This code was used to obtain the following results in the paper NMA and RW:

    a. The transition matrix for a RW on the aggregate network (Eq (26))
  
    b. The evidence flow network for comparison 1-3 of the Depression data shown in Figure 7
  
  
  C++ code that works out the flow of evidence for the comparison source-sink using the analytical RW approach
  
  The following data files must be provided by the user
  
  M.txt = list (length T) of the number of treatments in each trial
  
  SD.txt = list (length Q) of the standard deviation associated with each value of x
  
  B.txt = the edge-incidence matrix (QxN) of the network of treatments and trials, describes which treatments are in each trial
  
  Ba.txt = the edge-incidence matrix (KxN) of the aggregate network, describes which treatment comparisons are represented by each edge


  Using the data provided, a transition matrix for a RW on the aggregate network is derived
  
  Using the analogy to electrical networks, the code calculates the nodal voltages and edge currents for a circuit with a battery attached across source-sink where v(source)-1 and v(sink)=0
  
  The edge currents are then normalised by the total current leaving the source 
  
  These currents give the evidence flow from source to sink

Files:

Data files provided are for the Depression data described in Rucker and Schwarzer (2014). 

For this data:

N = 11 (number of nodes)

T = 26 (number of trials)

K = 20 (number of edges)

Q = 47 (total number of pairwise comparisons)

The files are:

1. M.txt = list (length T) of the number of treatments in each trial
2. x.txt = list (length Q) of the relative treatment effects measured in each trial
3. SD.txt = list (length Q) of the standard deviation associated with each value of x
4. B.txt = the edge-incidence matrix (QxN) of the network of treatments and trials, describes which treatments are in each trial
5. Ba.txt = the edge-incidence matrix (KxN) of the aggregate network, describes which treatment comparisons are represented by each edge

The trials are ordered in each file starting with 4-arm trials, then 3-arm trials and then 2-arm trials.

Using the trial labels in Rucker and Schwarzer, the order of trials in the files provided is:

58    (4-arm)

2     (3-arm)

16

24

43

49

82

86

104

5     (2-arm)

6

12

59

78

79

84

103

107

110

115

117

118

124

125

126

127
