//Code that implements the algorithm in Papakonstaninou et al (2018). 
//This version averages over the 'Random' version of the algorithm. 
//It is referred to as 'Average' in the NMA and RW manuscript. 
//The code relies on the Eigen library, the download files and documentation can be found here: http://eigen.tuxfamily.org/index.php?title=Main_Page

#include <iostream>
#include<Eigen/Dense>
#include<Eigen/LU>
#include<random>
#include<cmath>
#include<fstream>
#include<sstream>
#include<string> 
#include <algorithm>
#include <cstdlib>
#include<vector>
#include<iterator>
#include "time.h" 
using namespace std;
using namespace Eigen;


const int N{ 11 };//Number of nodes (treatments)
const int source = 1; //starting node
const int sink = 3; //end node
vector<vector<int>> all_paths; //vector of vectors to store all the path vectors found

//NB: nodes (treatments) are labelled as 1...N

//Recursive algorithm (Depth first search) to find all possible paths in a directed network
//algorithm based on video from Ashish Kumar April 2020: https://www.youtube.com/watch?v=TrEloQBv7WQ
void dfs(int visited[N], int u, int sink, vector<int> path, vector<vector<int>> poss_nodes) {
	visited[u-1] = 1; //u has been visited
	if (u == sink) {
		//once path has reached the sink - add path to list of paths
		all_paths.push_back(path);
		
	}
	else {
		//the path has not yet reached the sink
		//look for neighbours of u
		//number of neighbours = poss_nodes[u].size() 
		for (int i = 0; i < poss_nodes[u-1].size(); i++) {		
			int y = poss_nodes[u-1][i]; //y is the i'th neighbour of u 
			if (visited[y-1] == 0) { //if y has not yet been visited
				visited[y-1] = 1; //y has now been visited
				path.push_back(y); //add y to path
				dfs(visited, y, sink, path, poss_nodes); //keep looking for nodes until we reach the sink
				//now we look for other paths 
				path.pop_back();
			}
		}

	}

	visited[u-1] = 0; //find other paths from u
}




int main()
{
	time_t start, end; 
	time(&start);
	
	//read in the row of the h matrix that relates to the comparison source-sink
	//this should include all possible pairwise compairsons
	//if a certain pairwise compairson has no direct evidence set this H element = 0
	//comparison should be ordered as such: 12, 13, 14..1N, 23, 24,...2N...(N-1)N
	VectorXd h(N*(N-1)/2);
	int i0 = 0;
	double value0;
	ifstream file0;
	file0.open("Hflow13.txt");
	while (file0 >> value0)
	{
		h[i0] = value0;   //assigns to the array
		i0++;
	}

	cout << "h flow for 13" << endl << h << endl;

	
	//***************************************************************************************
	//Create transition matrix:**************************************************************
	//Create an NxN matrix where the diagonals=0 and element Q(i,j) = H_ij
	MatrixXd Q = MatrixXd::Zero(N, N); //initialise all equal to zero
	for (int i = 0; i < N; i++) {
		int pos = 0;
		for (int k = 0; k < i;k++) {
			pos = pos + N - (k + 1);
		}
		for (int j = 0; j < N; j++) {
			if (i != j) {
				if (i < j) {
					int index = pos + (j - i - 1);
					Q(i,j) = h[index];
				}
			}
		}
	}
	//H_ij = -H_ji 
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (i != j) {
				Q(j, i) = -Q(i, j);
			}
		}
	}


	//RW can only move in direction of flow
	//If H_ij<0 then set = 0
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (i != j) {
				if (Q(i, j) <= 0) {
					Q(i, j) = 0.0;
				}
				
			}
		}
	}

	//Create the transiton matrix by normalising the flow values in each row
	MatrixXd P = MatrixXd::Zero(N, N); //initialise all equal to zero
	for (int i = 0; i < N; i++) {
		double sum_row = 0.0;
		for (int j = 0; j < N; j++) {
			sum_row = sum_row + Q(i, j);
		}
		
		if (sum_row != 0) {
			for (int j = 0; j < N; j++) {
				P(i, j) = Q(i, j) / sum_row;
			}
		}
	}
	


	//Find all possible paths:****************************************************************************
	vector<vector<int>> poss_nodes; //vector of vectors storing the neighbours of each node
	
	for (int i = 0; i < N; i++) {		
		vector<int> poss_nodes_i;
		for (int j = 0; j < N; j++) {
			if (P(i, j) > 0) {
				//nodes are labelled 1,..,N (but j goes from 0...N-1)
				//hence we push back j+1
				poss_nodes_i.push_back(j + 1);				
			}
		}
		poss_nodes.push_back(poss_nodes_i);
	}
	
	

	
	int visited[N] = { 0 };//vector to store visited nodes	
	vector<int> path;//vector storing the nodes in a path
	path.push_back(source); //first node in path is the source

	//find all possible paths using a depth first search algorithm
	//this populates the vector of paths 'all_paths' 
	dfs(visited, source, sink, path, poss_nodes);

		
	//print path lengths and paths
	for (int i = 0; i < all_paths.size(); i++) {
		cout << "path " << i << ", length = " << all_paths[i].size()-1 << endl;
		for (int j = 0; j < all_paths[i].size(); j++) {
			cout << all_paths[i][j] << ", ";
		}
		cout << endl;
	}
	cout << endl;


	//Papakonstaninou algorithm 'Average':****************************************************************
	vector<double> sumProp(N*(N-1)/2,0.0);		//vector for working out average over proporiton contributions
	vector<double> sumPhi(all_paths.size(),0.0);	//vector for working out average over path flow

	

	int MAX = 100000; //number of iterations of Average
	for(int w=0; w<MAX; w++){
		cout << "iteration = " << w << endl;
		
		//at each run we must reset the graph and its matrix of flow***********************
		MatrixXd Q = MatrixXd::Zero(N, N); //initialise all equal to zero
		for (int i = 0; i < N; i++) {
			int pos = 0;
			for (int k = 0; k < i;k++) {
				pos = pos + N - (k + 1);
			}
			for (int j = 0; j < N; j++) {
				if (i != j) {
					if (i < j) {
						int index = pos + (j - i - 1);
						Q(i,j) = h[index];
					}
				}
			}
		}

		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				if (i != j) {
					Q(j, i) = -Q(i, j);
				}
			}
		}


		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				if (i != j) {
					if (Q(i, j) <= 0) {
						Q(i, j) = 0.0;
					}
					
				}
			}
		}
		//***********************************************************************
		

		//choose a random path:
			//pick an i from 0 to all_paths.size()-1
			//then chosen path is = all_paths[i]
			//length of this path (no. of edges) = all_paths[i].size()-1 
			//for j=0->all_paths[i].size()-1: all_paths[i][j] is the set of nodes

		double sumQ = 100; //sumQ keeps track of whether all the edge flows has been assigned to paths
		//random device for path selection
		random_device gen_u;
		mt19937 twist_u(gen_u());
		
		//keep track of which paths have been selected (reset vectors)
		vector<int> RandPath;  //vector storing the labels (1...no. of paths) of paths that have been selected
		vector<double> RandFlow; //vector storing the flow assigned to each selected path
		

		//while all edge flow has not been assigned to path flow
		while(sumQ>1e-10){ //set small value to deal with precision of computer

			int i;
			int check=0;
			while(check==0){ //while the path has already been selected...
				check = 1;
				//select a path from all_paths at random
				uniform_int_distribution<int> uni_dist_u(0, all_paths.size()-1);
				i = uni_dist_u(twist_u);
				for(int j=0; j<RandPath.size(); j++){ //check if path has already been selected					
					if(RandPath[j]==i){ check = 0; } //this path has already been selected
				}
				
			}
			RandPath.push_back(i); //add path to selected paths
			int li = (int)all_paths[i].size() - 1; //length of randomly selected path (no. of edges = no. of nodes - 1)
			
			
			//flow along each edge in the path - from matrix of flow
			vector<double> Fi(li);
			for(int k=0; k<li; k++){
				Fi[k] = Q(all_paths[i][k]-1,all_paths[i][k+1]-1);
			}

			
			//find min flow through edges 
			double F_min = *min_element(Fi.begin(), Fi.end());
			RandFlow.push_back(F_min); //flow assigned to path i

			//subtract path flow from each edge in path i
			for(int k=0; k<li; k++){
				Q(all_paths[i][k]-1,all_paths[i][k+1]-1) = Q(all_paths[i][k]-1,all_paths[i][k+1]-1) - F_min;
				
			}
			
			//New matrix of flow Q:
			//work out total flow remaining in the network
			sumQ = 0.0;
			for (int p = 0; p < N; p++) {
				for (int q = 0; q < N; q++) {
					
					sumQ = sumQ + Q(p, q);	
				}
			}
			
			//if there is still flow left then select another path

		}


		

		//Streams identified in iteration w************************************************
		for(int j=0; j<RandPath.size(); j++){
			int path_label = RandPath[j];	
			sumPhi[path_label] = sumPhi[path_label]+RandFlow[j];	//sum flow at each iteration to obtain average
			/*cout << "Path number = " << RandPath[j] << ", flow = " << RandFlow[j] << endl << "Path: ";
			for (int k = 0; k < all_paths[path_label].size(); k++) {
				cout << all_paths[path_label][k] << ", ";
			}
			cout << endl;*/
		}
		


		//Find the contribution of each edge***************************************************************************************
		vector<double> prob_edge(N*(N-1)/2, 0.0);
		int ind = 0; //index labelling the pairwise comparison (from 12, 13,...,(N-1)N)
		for (int i = 1; i < N+1; i++) {
			for (int j = i + 1; j < N+1; j++) {
				//edge ij
				for (int k = 0; k < RandPath.size(); k++) {
					int path_label = RandPath[k];
					//does path path_label contain the edge ij?
					int edgeij = 0;
					for (int l = 0; l < all_paths[path_label].size()-1; l++) {
						if (all_paths[path_label][l] == i && all_paths[path_label][l+1] == j) {
							edgeij = 1; //yes it does contain edge ij
						}
						else if (all_paths[path_label][l] == j && all_paths[path_label][l + 1] == i) {
							edgeij = 1; //yes it does contain edge ij (ji)
						}
					}
					//if yes: p_ab = sum_{paths i containing ab} phi_i/length(pi_i)
					if (edgeij == 1) {
						double length = (double)(all_paths[path_label].size()-1); 
						//add phi_i/length(pi_i) to p_ab
						prob_edge[ind] = prob_edge[ind] + (RandFlow[k] / length);
					}
				}
				
				sumProp[ind] = sumProp[ind]+prob_edge[ind];//add to vector to work out average over iterations
				ind = ind + 1;
			}
			
		}

	} //end of iterations*******************************************************************

	//work out average prop contribution*********************************************************************

	vector<double> MeanProp(N*(N-1)/2);
	cout << "average over prop contributions" << endl;
	for(int i=0; i<N*(N-1)/2; i++){
		MeanProp[i] = sumProp[i]/((double)MAX);
		cout << MeanProp[i] << endl;
		
	}

	//work out average path flow in each possible path*********************************************************************
	vector<double> MeanPhi(all_paths.size());
	cout << "average over phi" << endl;
	for(int i=0; i<all_paths.size(); i++){
		MeanPhi[i] = sumPhi[i]/((double)MAX);
		cout << "Path label = " << i << ", Mean phi = " << MeanPhi[i] << endl;
		
	}

	//check average prop contribution again prop contribution found from average path flow*************************************
	//(these should give the same result)
	vector<double> AVprob_edge(N*(N-1)/2, 0.0);
	cout << "resulting average over prop: " << endl;
	int ind = 0;
	for (int i = 1; i < N+1; i++) {
		for (int j = i + 1; j < N+1; j++) {
			cout << "edge: " << i << j;
			for (int k = 0; k < all_paths.size(); k++) {
				//does path k contain the edge ij?
				int edgeij = 0;
				
				for (int l = 0; l < all_paths[k].size()-1; l++) {
					if (all_paths[k][l] == i && all_paths[k][l+1] == j) {
						edgeij = 1;//yes it does contain edge ij
					}
					else if (all_paths[k][l] == j && all_paths[k][l + 1] == i) {
						edgeij = 1; //yes it does contain edge ij (ji)
					}
				}
				//yes:
				if (edgeij == 1) {
					double length = (double)(all_paths[k].size()-1);
					AVprob_edge[ind] = AVprob_edge[ind] + MeanPhi[k] / length;
				}
			}
			cout << ": contribution = " << AVprob_edge[ind] << endl;
			ind = ind + 1;
		}
		
	}

	

	time(&end);

	double time_taken = double(end-start); 
    	cout << "Time taken by program is : " << fixed << time_taken; 
    	cout << " sec " << endl; 

	return 0;

}
