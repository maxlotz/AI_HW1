/*Compiler settings
For C++, Kattis uses gcc version g++ (Ubuntu 5.4.0-6ubuntu1~16.04.2) 5.4.0 20160609 with the following flags: -g -O2 -static -std=gnu++11 {files}.
System libraries
You are allowed to use all standard libraries included with C++.
Test with ./a.out <samples/hmm3_01.in */

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <cmath>

using namespace std;

typedef vector< vector<double> > matrix; // matrices are all vectors of vectors of doubles

struct s1
{
	vector<matrix> mat3d; // essentially a 3d array, its a vector<vector<vector<double>>>
	matrix mat1;
	matrix mat2;
	matrix mat3;
	vector<double> vec;
};

// FUNCTION DECLARATIONS AND EXPLANATIONS  HERE

// Displays the contents of a matrix to the terminal, columns seperated by tab, rows seperated by newline
void dispmat(matrix inmat); 

// Returns transpose of given matrix
matrix transpose (matrix inmat);

// Converts string from input file into a matrix
matrix str2mat (string line);

// Converts string from input file into vector of ints containing the emmission sequence
vector<int> str2seq (string line);

// Converts matrix into string for output format (opposite of str2mat)
string mat2str (matrix inmat);

// Returns probability of getting first observation of the emmission sequence
matrix alpha1 (matrix B, matrix PI, vector<int> O_seq);

// Returns struct containing a matrix (row vector) of alpha values (probability of given observation occuring at that timestep) and a matrix (row vector) of their maximum index (most likely state at that timestep)
s1 forward_pass(matrix A, matrix B, matrix PI, vector<int> O_seq);

// Performs backward pass and returns beta matrix
matrix backward_pass(vector<double> c, matrix A, matrix B, vector<int> O_seq);

// performs gamma pass, calculates gamma_i (vector) and gamma_ij (matrix) and returns as object of s1
s1 gammas (matrix alphas, matrix betas, matrix A, matrix B, vector<int> O_seq);

// Performs re-estimation
s1 Re_estimate(vector<matrix> gamma_ij, matrix gamma_i, vector<int> O_seq);

// computelogprob
double computelogprob(vector<double> c);

// Made timesteps as global variables here.
int T;
int N;
int K;
int maxiters = 100;

// MAIN PROGRAM HERE
int main(void)
{
  	string line;
	matrix A; //transmission matrix
	matrix B; //emmission matrix
	matrix PI; //intial state distribution
	vector<int> O_seq; //emmission sequence
   	
	getline(cin,line);	
	A = str2mat(line);
	getline(cin, line);
	B = str2mat(line);
	getline(cin, line);
	PI = str2mat(line);
	getline(cin, line);
	O_seq = str2seq(line);
	
	T = O_seq.size();		// no. of timesteps (observation sequence size)
	N = A.size();			// no. of states
	K = B[0].size();		// no. of observations
	
	s1 out;
	s1 out2;
	s1 out3;
	matrix beta;

	out = forward_pass (A, B, PI, O_seq);
	beta = backward_pass(out.vec, A, B, O_seq);
	out2 = gammas (out.mat1, beta, A, B, O_seq);
	out3 = Re_estimate(out2.mat3d, out2.mat1, O_seq);

	cout << "\n\nA =\n\n";
	dispmat(out3.mat1);	

	cout << "\n\nB =\n\n";
	dispmat(out3.mat2);	

	cout << "\n\nPI =\n\n";
	dispmat(out3.mat3);	
	return 0; 
}

// FUNCTIONS BODIES HERE
void dispmat(matrix inmat)
{	
	int a = T + 2;
	for (int i = 0; i < inmat.size(); i++)
	{
	    for (int j = 0; j < inmat[i].size(); j++)
	    {
	        cout << inmat[i][j] << "\t";
	    }
	    cout << "\n";
	}
}

matrix transpose (matrix inmat)
{
	matrix outmat (inmat[0].size(), vector<double> (inmat.size()) );
	for (int i = 0; i < inmat.size(); i++)
	{
		for (int j = 0; j < inmat[0].size(); j++)
		{
			outmat[j][i] = inmat[i][j];
		}
	}
	return outmat;
}

matrix str2mat (string line)
{
	istringstream ss(line);
	vector <string> record;

	while(ss)
	{
		string s;
		{
			if (!getline( ss, s, ' ' )) break;
  			record.push_back( s );
		}
	}
	int rows = stoi(record[0]);
	int cols = stoi(record[1]);
	matrix outmat (rows, vector<double> (cols));

	for(int i = 0; i < rows; i++)
	{
		for(int j = 0; j < cols; j++)
		{
			outmat[i][j] = stod(record[i*cols + j + 2]);
		}
	}
	return outmat;
}

vector<int> str2seq (string line)
{
	istringstream ss(line);
	vector <string> record;

	while(ss)
	{
		string s;
		{
			if (!getline( ss, s, ' ' )) break;
  			record.push_back( s );
		}
	}
	int vecsize = stoi(record[0]);
	vector<int> outvec(vecsize);
	
	for (int i = 0; i < vecsize; i++)
	{
		outvec[i] = stoi(record[i+1]);
	}
	return outvec;
}

string mat2str (matrix inmat)
{
	int rows = inmat.size();
	int cols = inmat[0].size();
	string outline;
	string blank = " ";
	outline += to_string(rows) + blank + to_string(cols);

	for(int i = 0; i < rows; i++)
	{
		for(int j = 0; j < cols; j++)
		{
			stringstream ss;
			ss << inmat[i][j];
			string temp = ss.str();
			if (temp == "0")
			{
				temp = "0.0";
			}
			outline += blank + temp;
		}
	}
	return outline;
}

s1 forward_pass (matrix A, matrix B, matrix PI, vector<int> O_seq)
{
	s1 out;
	// creates alpha and c matrix of correct size and fills with zeros
	matrix alpha (T, vector<double>(N));
	vector<double> c(T);

	// compute alpha0(i)
	for (int i = 0; i < N; i++)	// from 0 to N-1
	{
		alpha[0][i] = PI[0][i]*B[i][O_seq[0]];
		c[0] += alpha[0][i];
	}
	// scale alpha0(i)
	c[0] = 1/c[0];
	for (int i = 0; i < N; i ++)
	{
		alpha[0][i] = c[0]*alpha[0][i];
	}
	// compute alphat(i)
	for (int t = 1; t < T; t++) // from 1 to T-1
	{
		for (int i = 0; i < N; i++) // from 0 to N-1
		{
			for (int j = 0; j < N; j++) // from 0 to N-1
			{
				alpha[t][i] += alpha[t-1][j]*A[j][i];
			}
			alpha[t][i] = alpha[t][i]*B[i][O_seq[t]];
			c[t] += alpha[t][i];
		}
		// scale alphat(i)
		c[t] = 1/c[t];
		for (int i = 0; i < N; i++)
		{
			alpha[t][i] = c[t]*alpha[t][i];
		}
	}
	out.mat1 = alpha;
	out.vec = c;
	return out;
}

matrix backward_pass(vector<double> c, matrix A, matrix B, vector<int> O_seq)
{
	matrix beta (T, vector<double>(N)); // creates beta matrix size TxN filled with zeros
	// beta_t-1(i) = 1 scaled by C_t-1
	for (int i = 0; i < N; i ++) beta[T-1][i] = c[T-1];
	// backward pass
	for (int t = T-2; t > -1; t--) // from T-2 to 0
	{
		for (int i = 0; i < N; i ++) // from 0 to N-1
		{
			for (int j = 0; j < N; j++) // from 0 to N-1
			{
				beta[t][i] += A[i][j]*B[j][O_seq[t+1]]*beta[t+1][j];
			}
			beta[t][i] = c[t]*beta[t][i];
		}
	}
	return beta;
}

s1 gammas (matrix alphas, matrix betas, matrix A, matrix B, vector<int> O_seq)
{
	s1 out;
	matrix gamma_i (T, vector<double>(N)); // initialises empty matrix size TxN with zeros.
	vector<matrix> gamma_ij(T-1, matrix(N, vector<double>(N))); // initialises empty 3d matrix size T-1xNxN with zeros.

	// calculates gamma_ij amd gamma_i for t=0 to T-2, making them both length T-1
	for (int t = 0; t < T-1; t++) // from 0 to T-1
	{
		double denom = 0;
		for (int i = 0; i < N; i ++) // from 0 to N-1
		{
			for(int j = 0; j < N; j++)
			{
				denom += alphas[t][i]*A[i][j]*B[j][O_seq[t+1]]*betas[t+1][j];
			}
		}
		for (int i = 0; i < N; i++) // from 0 to N-1
		{
			for (int j = 0; j < N; j++)
			{
				gamma_ij[t][i][j] = (alphas[t][i]*A[i][j]*B[j][O_seq[t+1]])*betas[t+1][j]/denom;
				gamma_i[t][i] += gamma_ij[t][i][j]; 
			}
		}
	}
	// special case, final value for gamma_i, now it is size T
	// This definitely works correctly
	double denom = 0;
	for (int i = 0; i < N; i++)
	{
		denom += alphas[T-1][i];
	}
	for (int i = 0; i < N; i++)
	{
		gamma_i[T-1][i] = alphas[T-1][i]/denom;
	}

	out.mat3d = gamma_ij;
	out.mat1 = gamma_i;
	return out;
}

s1 Re_estimate(vector<matrix> gamma_ij, matrix gamma_i, vector<int> O_seq)
{	
	s1 out;
	// Initialises empty matrices
	matrix A (N, vector<double>(N));
	matrix B (N, vector<double>(K));
	matrix PI (1, vector<double>(N));

	//re-estimate PI
	for(int i = 0; i < N; i++) PI[0][i] = gamma_i[0][i]; // from O to N-1

	//re-estimate A
	for(int i = 0; i < N; i++)	// from 0 to N-1
	{
		for(int j = 0; j < N; j++) // from 0 to N-1
		{
			double numer, denom = 0;
			for(int t = 0; t < T-1; t++) // from 0 to T-2
			{
				numer+= gamma_ij[t][i][j];
				denom+= gamma_i[t][i];
			}
			A[i][j] = numer/denom;
		}
	}

	// re-estimate B
	for (int i = 0; i < N; i++) // from 0 to N-1
	{
		for(int j = 0; j < N; j++) // from 0 to N-1
		{
			double numer, denom = 0;
			for(int t = 0; t < T; t++) // from 0 to T-1
			{
				if(O_seq[t] == j) numer+= gamma_i[t][i];
				denom += gamma_i[t][i];
			}
			B[i][j] = numer/denom;
		}
	}

	out.mat1 = A;
	out.mat2 = B;
	out.mat3 = PI;
	return out;
}

double computelogprob(vector<double> c)
{
	double logprob = 0;
	for(int i = 0; i < T; i++) logprob+= log(c[i]);
	return -logprob;
}


