//Code to perform the graph theoretical model NMA
//A one-step and two-step version are provided (the one step version is commented out)
//produces the hat matrix of the aggregate model (used to obtain evidence flow networks)
//values defined are for the Depression data set in Rucker and Schwarzer (2014)
//We perform a fixed effect model
//for a random effects model simply estimate tau (e.g. from method of moments)
//then add tau^2 to each element of the vector of variances in each trial Var and proceed as normal
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
#include<iomanip>

using namespace std;
using namespace Eigen;

const int N {11}; //number of treatments
const int T {26}; //number of trials
const int K {20}; //number of edges

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


int main(){
	//NB: my files are ordered differently than the order of trials given in Rucker and Schwarzer (2014)
	//They are ordered such that they start with 4-arm trials, then 3-arm, then 2-arm
	//If all files are defined consistently with each other this code works with the trials in any order

	
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

	//read in vector of treatment effects from each trial
	VectorXd x(Q);  
	int x0 = 0;
	double valuex0;
	ifstream filex0;
	filex0.open("x.txt");
	while (filex0 >> valuex0)
	{
		x[x0] = valuex0;   
		x0++;
	}
	//cout << "treatment effects" << endl << x << endl; 

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
	//***********************************************************************************************
	
	
	
	//Perform the ONE STEP model:*********************************************************************
	//Now work out hat matrix of the one-step model
/*	MatrixXd L(N,N);
	L = B.transpose()*W*B;
	cout << "Laplacian = " << endl << L << endl;
	MatrixXd L_pseudo(N,N);
	MatrixXd O = MatrixXd::Ones(N,N);
	MatrixXd A(N,N);
	A = L - O/(double)N;
	L_pseudo = A.inverse() + O/(double)N;
	cout << "Pseudo-inverse Laplacian = " << endl << L_pseudo << endl;

	MatrixXd H(Q,Q); 
	H = B*L_pseudo*B.transpose()*W;
	cout << "Hat matrix = " << endl << H << endl;

	//And now get the network treatment effect estimates from the one-step model
	VectorXd theta(Q);
	theta = H*x;
	cout << "theta = " << endl << theta << endl << endl;

	VectorXd expTheta(Q);
	for(int i=0; i<Q; i++){
		expTheta[i] = exp(theta[i]);
	}
	cout << "OR = " << endl << expTheta << endl;
*/

	//Perform the TWO STEP (aggregate model):***********************************************************
	cout << "TWO STEP:" << endl;
	
	VectorXd Xpair(K); //Vector of direct estimates
	MatrixXd Wa = MatrixXd::Zero(K,K); //matrix of aggregate weights

	//STEP ONE--------------------------------------------------------------------------------
	//perform a pairwise MA across each edge of the network using the adjusted edge weights
	//iterate over rows of Ba
	for(int i=0; i<K; i++){	
		int indexA; //element of row i of Ba that corresponds to baseline treatment
		int indexB; //element of row i of Ba that corresponds to the other treatment
		for(int j=0; j<N; j++){
			if(Ba(i,j)==1){ indexA = j; }
			else if(Ba(i,j)==-1){ indexB = j; }
		}
		//Now find all rows of B that have the same elements and perform a weighted mean
		double sumAB = 0.0;
		double WAB = 0.0;
		for(int k=0; k<Q; k++){
			if(B(k,indexA)==1 && B(k,indexB)==-1){
				//work out the weighted mean and inverse-variance of the weighted mean
				sumAB = sumAB + x[k]*W(k,k);
				WAB = WAB + W(k,k);
			}
		}
		Xpair[i] = sumAB/WAB;
		Wa(i,i) = WAB;

	}

	//STEP TWO--------------------------------------------------------------------------------
	//Now work out hat matrix of the aggregate model
	MatrixXd La(N,N);
	La = Ba.transpose()*Wa*Ba;
	//cout << "Aggregate laplacian La = " << endl << La << endl;
	MatrixXd La_pseudo(N,N);
	MatrixXd Oa = MatrixXd::Ones(N,N);
	MatrixXd Aa(N,N);
	Aa = La - Oa/(double)N;
	La_pseudo = Aa.inverse() + Oa/(double)N;
	//cout << "Aggregate Pseudo-inverse Laplacian = " << endl << La_pseudo << endl;

	MatrixXd Ha(K,K); 
	Ha = Ba*La_pseudo*Ba.transpose()*Wa;
	cout << "Aggregate Hat matrix = " << endl << Ha << endl;

	//Now work out treatment effect estimates from the aggregate model
	VectorXd Theta_a(K);
	Theta_a = Ha*Xpair;
	//cout << "network estimates from 2-step: " << endl << Theta_a << endl;

	VectorXd expTheta_a(K);
	for(int i=0; i<K; i++){
		expTheta_a[i] = exp(Theta_a[i]);
	}
	cout << "2step OR = " << endl << expTheta_a << endl;


	return 0;
}
