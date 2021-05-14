//Code to work out flow of evidence using the analytical RW approach
//read in data and calculate the aggregate weight matrix
//use this to define a transition matrix on the aggregate network
//from this derive the nodal potentials (with v_source=1 and v_sink=0)
//work out edge currents for v_source=1 then normalise currents to get edge currents when total current flowing out of source =1
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


const int N{ 11 }; //number of nodes (treatments)
const int T {26}; //number of trials
const int K {20}; //number of edges

const int source = 1; //starting node
const int sink = 3; //end node 
vector<vector<int>> all_paths;


//A function to work out the adjusted weights for a m-arm trial
//Var is a q=m(m-1)/2 length matrix containg the variances for each comparison in this trial
//m = no. of treatments in this trial and 'trial' is the trial label
//an m-arm trial describes a subgraph
MatrixXd adjust_weights(VectorXd Var, int m, int trial){
	
	//create variance matrix and B matrix for this subgraph 
	int q = (m*(m-1))/2; 
	MatrixXd B = MatrixXd::Zero(q,m);
	MatrixXd V = MatrixXd::Zero(m,m);
	int ind = 0;
	for(int k=1; k<m; k++){
		
		for(int i=0; i<m-k; i++){
			int index = ind + i;
			V(k-1,k+i) = Var[index];
			V(k+i,k-1) = Var[index];
			B(index,k-1) = 1.0;
			B(index,k+i) = -1.0;
		}
		
		ind = ind + (m-k);
		
	}	
	//cout << "V" << endl << V << endl;
	//cout << "B" << endl << B << endl;
	

	//calculate pseudo inverse of Laplacian for this subgraph
	MatrixXd L_pseudo(m,m);
	L_pseudo = -(1.0/(2.0*pow((double)m, 2.0)))*B.transpose()*B*V*B.transpose()*B;
	//cout << "L pseudo inverse " << endl << L_pseudo << endl;		
	MatrixXd O = MatrixXd::Ones(m,m);
	//calculate the laplacian	
	MatrixXd L(m,m);
	MatrixXd LO(m,m);
	LO = L_pseudo - (1.0/(double)m)*O;
	L = LO.inverse() + (1.0/(double)m)*O;	
	//cout << "L" << endl << L << endl;
	
	//create weight matrix (for this subgraph) from L
	MatrixXd W = MatrixXd::Zero(q,q);
	ind = 0; //reset ind
	for(int k=1; k<m; k++){
		for(int i=0; i<m-k; i++){
			int index = ind + i;
			W(index,index)=-L(k-1,k+i);
		}
		ind = ind + (m-k);
	}
	//cout << "W" << endl << W << endl;

	return W;

}

//Function to remove a row from a matrix
void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove)
{
	unsigned int numRows = matrix.rows()-1;
	unsigned int numCols = matrix.cols();

	if( rowToRemove < numRows ){
		matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.bottomRows(numRows-rowToRemove);
	}


	matrix.conservativeResize(numRows,numCols);
}

//Function to remove a column from a matrix
void removeColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove)
{
	unsigned int numRows = matrix.rows();
	unsigned int numCols = matrix.cols()-1;

	if( colToRemove < numCols ){
		matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.rightCols(numCols-colToRemove);
	}

	matrix.conservativeResize(numRows,numCols);
}

int main()
{
	//Read in Data*************************************************************************************
	//read in vector M (number of treatments in each trial)
	VectorXi m(T); 
	int m0 = 0;
	int valuem0;
	ifstream filem0;
	filem0.open("M.txt");
	while (filem0 >> valuem0)
	{
		m[m0] = valuem0;   
		m0++;
	}

	//create vector q = number of pairwise comparisons in each trial
	VectorXi q(T); 
	int Q = 0; //number of pairwise comparisons
	for(int i=0; i<T; i++){
		q[i] = (m[i]*(m[i]-1))/2;
		Q = Q + q[i];
	}
	//cout << "m" << endl << m << endl << "q" << endl << q << endl << "Q = " << Q << endl;

	//read in vector of variances from each trial
	VectorXd Var(Q);
	int v0 = 0;
	double valuev0;
	ifstream filev0;
	filev0.open("SD.txt");
	while (filev0 >> valuev0)
	{
		Var[v0] = pow(valuev0,2.0);   //square SD to get variance
		v0++;
	}
	//cout << "variances " << endl << Var << endl;	

	//read in B matrix (describes which treatments are in each trial)
	//B has Q rows and N columns
	//file should be structured as:
	//B_00	B_01	B_02...
	//B_10	B_11	B_12...
	//...
	//columns should be in order 12, 13, 14, 23, 24... for each trial
	ifstream fileB("B.txt"); 
	MatrixXd B(Q,N);
	for(int i=0; i<Q; i++){
		for(int j=0; j<N; j++){
			fileB >> B(i,j);
		}
	}
	//cout << "B" << endl << B << endl;
	
	//read in aggregate Ba (describes which treatments are compared in each edge)
	//Ba has K rows (no. of edges) and N columns (no. of nodes)
	//file should be structured as:
	//Ba_00	Ba_01	Ba_02...
	//Ba_10	Ba_11	Ba_12...
	//...
	//columns should be in order 12, 13, 14, 23, 24... 
	ifstream fileBa("Ba.txt"); 
	MatrixXd Ba(K,N);
	for(int i=0; i<K; i++){
		for(int j=0; j<N; j++){
			fileBa >> Ba(i,j);
		}
	}
	//cout << "Ba" << endl << Ba << endl;
	//Finished reading in data*********************************************************************
	

	//adjust weights to deal with correlations from multi-arm trials ******************************
	MatrixXd W = MatrixXd::Zero(Q,Q); //full weight matrix
	int sum_q=0;
	for(int i=0; i<T; i++){
		//cout << "trial " << i+1 << endl;
		VectorXd Vari = VectorXd::Zero(q[i]);
		
		for(int j=0; j<q[i]; j++){
			Vari[j] = Var[sum_q+j];
		}
		MatrixXd Wi = MatrixXd::Zero(q[i],q[i]);
		Wi = adjust_weights(Vari, m[i], i+1);
		//cout << "Wi = " << endl << Wi << endl;

		for(int k=0; k<q[i]; k++){
			W(sum_q+k,sum_q+k) = Wi(k,k);
		}
		sum_q = sum_q + q[i];

	}
	//cout << "full (adjusted) weight matrix = " << endl << W << endl;
	
	
	//Work out aggregate weight matrix***************************************************************
	MatrixXd Wa = MatrixXd::Zero(K,K); //matrix of aggregate weights

	//calculate inverse variance of weighted mean
	//iterate over rows of Ba
	for(int i=0; i<K; i++){	
		int indexA; //element of row i of Ba that corresponds to baseline treatment
		int indexB; //element of row i of Ba that corresponds to the other treatment
		for(int j=0; j<N; j++){
			if(Ba(i,j)==1){ indexA = j; }
			else if(Ba(i,j)==-1){ indexB = j; }
		}
		//Find all rows of B that have the same elements and calculate inverse variance of weighted mean
		double WAB = 0.0;
		for(int k=0; k<Q; k++){
			if(B(k,indexA)==1 && B(k,indexB)==-1){
				WAB = WAB + W(k,k);
			}
		}
		Wa(i,i) = WAB;
	}

	

	//create transition matrix***************************************************************
	//using aggregate weights define a transition matrix for a walker on the aggregate network
	MatrixXd Q1 = MatrixXd::Zero(N, N); //initialise all equal to zero
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			double Wij=0.0;
			//search Ba:
			for(int k=0; k<K; k++){
				if(Ba(k,i)==1 && Ba(k,j)==-1){
					Wij = Wa(k,k);
				}
				else if(Ba(k,i)==-1 && Ba(k,j)==1){
					Wij = Wa(k,k);
				}
			}
			Q1(i,j) = Wij;
		}
	}
	//cout << "Q1" << endl << Q1 << endl;

	//normalise each row of Q1 to obtain transition matrix
	MatrixXd P = MatrixXd::Zero(N, N); //initialise all equal to zero
	for (int i = 0; i < N; i++) {
		double sum_row = 0.0;
		for (int j = 0; j < N; j++) {
			sum_row = sum_row + Q1(i, j);
		}
		if (sum_row != 0) {
			for (int j = 0; j < N; j++) {
				P(i, j) = Q1(i, j) / sum_row;
			}
		}
	}
	//edit row corresponding to sink:
	for(int j=0; j<N; j++){
		if(j==sink-1){
			P(sink-1,j) = 1.0;
		}
		else{
			P(sink-1,j) = 0.0;
		}
	}	
	//cout << "P:" << endl << P << endl;

	MatrixXd Pfull = MatrixXd::Zero(N, N); //full transition matrix to hold on to
	for (int i = 0; i < N; i++) {		
		for (int j = 0; j < N; j++) {
			Pfull(i, j) = P(i,j);
		}
		
	}

	//now edit P for reduced transition matrix***************************************************************
	//remove rows and columns referring to source and sink
	//remove the node correpsonding to the later column/row first

	//remove rows
	if(sink>source){
		removeRow(P,sink-1);
		removeRow(P,source-1);
	}
	else{
		removeRow(P,source-1);
		removeRow(P,sink-1);
	}
	//get vector of transitions from source (column correpsonding to source)	
	VectorXd P_1 = VectorXd::Zero(N - 2); //column of P referring to inflows to source (1) with source and sink row removed
	for(int i=0; i<N-2; i++){
		P_1[i] = P(i,source-1);
	}

	//remove columns	
	if(sink>source){
		removeColumn(P,sink-1);
		removeColumn(P,source-1);
	}
	else{
		removeColumn(P,source-1);
		removeColumn(P,sink-1);
	}

	cout << "P full " << endl << Pfull << endl;
	cout << "P reduced" << endl << P << endl;
	cout << "P1 " << endl << P_1 << endl;
	
	MatrixXd I = MatrixXd::Identity(N - 2, N - 2); //N-2 x N-2 identity matrix

	//Calculate nodal potentials from harmonic function********************************************************
	VectorXd v_red(N - 2); //vector of unknown potentials
	v_red = (I - P).inverse()*P_1;
	cout << "nodal voltages (reduced):" << endl << v_red << endl;

	
	//now get full vector of potentials (v[sink-1]=0 v[source-1]=1)
	VectorXd v(N);
	for(int i=0; i<N; i++){
		if(i<source-1 && i<sink-1){
			v[i] = v_red[i];
		}
		else if(source<sink){
			if(i==source-1){
				v[i] = 1.0;
			}
			else if(i>source-1 && i<sink-1){
				v[i] = v_red[i-1];
			}
			else if(i==sink-1){
				v[i] = 0.0;
			}
			else if(i>source-1 && i>sink-1){
				v[i] = v_red[i-2];
			}
			else{
				cout << "i out of range" << endl;
			}
		}
		else{ //sink<source
			if(i==sink-1){
				v[i] = 0.0;
			}
			else if(i>sink-1 && i<source-1){
				v[i] = v_red[i-1];
			}
			else if(i==source-1){
				v[i] = 1.0;
			}
			else if(i>sink-1 && i>source-1){
				v[i] = v_red[i-2];
			}
			else{
				cout << "i out of range" << endl;
			}
		}
	}
	cout << "full v " << endl << v << endl;

	//work out edge currents for v_source = 1********************************************************
	VectorXd I_v1(K);
	I_v1 = Wa*Ba*v;
	cout << "currents at v=1 " << endl << I_v1 << endl;

	//Total current out of source
	//column of Ba correpsonding to source = Ba(i,source-1)
	double I_tot = 0.0;
	for(int i=0; i<K; i++){
		if(Ba(i,source-1)!=0){
			I_tot = I_tot + abs(I_v1[i]); //current out of source is positive
		}
	}

	//find edge currents when current out of source = 1 (i.e. normalise I_v1)************************
	VectorXd I_i1(K); //normalised currents = evidence flow
	for (int k = 0; k < K; k++) {
		I_i1[k] = I_v1[k] / I_tot;
	}
	//these are the evidence flows for comparison source-sink
	cout << "currents with Ia=1:" << endl << I_i1 << endl; 



	return 0;

}