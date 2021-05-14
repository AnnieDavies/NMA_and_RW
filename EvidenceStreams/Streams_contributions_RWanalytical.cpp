//Code to work out streams of evidence and proportion contributions using the analytical RW approach
//Define a transition matrix on the evidence flow network using a row of the H matrix
//find all possible paths on the evidence flow network 
//flow in path = product over transition probabilities in each edge along the path
//use these to calculate proportion contributions
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

using namespace std;
using namespace Eigen;


const int N{ 11 };//number of nodes (treatments)
const int source = 1; //starting node
const int sink = 3; //end node
vector<vector<int>> all_paths; //vector of vectors to store all the path vectors found


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
	//User defined variables:***************************************************************
	
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
	
	cout << "P:" << endl << P << endl;

	//Find all possible paths:****************************************************************************
	vector<vector<int>> poss_nodes; //vector of vectors storing neighbours of each node
	
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


	//now work out the probability of each path***************************************************************************
	//P(path) = product over transition probability in each edge
	vector<double> prob_paths;
	for (int i = 0; i < all_paths.size(); i++) {
		double prob = 1.0;
		for (int j = 0; j < all_paths[i].size(); j++) {
			cout << all_paths[i][j] << ", ";
		}
		for (int j = 0; j < all_paths[i].size()-1; j++) {
			prob = prob * P(all_paths[i][j]-1, all_paths[i][j + 1]-1);
		}
		cout << "  prob = " << prob << endl;
		prob_paths.push_back(prob);
		
	}

	
	//Find the contribution of each edge***************************************************************************************
	vector<double> prob_edge(N*(N-1)/2, 0.0);
	int ind = 0;
	for (int i = 1; i < N+1; i++) {
		for (int j = i + 1; j < N+1; j++) {
			//edgeij
			cout << "edge: " << i << j;
			for (int k = 0; k < all_paths.size(); k++) {
				//does path k contain the edge ij?
				int edgeij = 0;
				for (int l = 0; l < all_paths[k].size()-1; l++) {
					if (all_paths[k][l] == i && all_paths[k][l+1] == j) {
						edgeij = 1; //yes it does contain edge ij
					}
					else if (all_paths[k][l] == j && all_paths[k][l + 1] == i) {
						edgeij = 1; //yes it does contain edge ij (ji)
					}
				}
				//if yes: p_ab = sum_{paths i containing ab} phi_i/length(pi_i)
				if (edgeij == 1) {
					double length = (double)(all_paths[k].size()-1);
					prob_edge[ind] = prob_edge[ind] + prob_paths[k] / length; //add phi_i/length(pi_i) to p_ab
				}
			}
			cout << ": contribution = " << prob_edge[ind] << endl;
			ind = ind + 1;
		}
		
	}

	

	return 0;

}
