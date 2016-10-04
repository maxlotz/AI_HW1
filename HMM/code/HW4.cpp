/*Compiler settings

For C++, Kattis uses gcc version g++ (Ubuntu 5.4.0-6ubuntu1~16.04.2) 5.4.0 20160609 with the following flags: -g -O2 -static -std=gnu++11 {files}.
System libraries

You are allowed to use all standard libraries included with C++.

Test with ./a.out <samples/hmm4_01.in */

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <cmath>
#include <float.h>

using namespace std;

typedef vector< vector<double> > matrix; // matrices are all vectors of vectors of doubles

class HMM {
	private:

	public:
		matrix A, B, PI, alpha, beta, gamma_i;
		vector<int> O_seq;
		vector<double> c, storeprob;
		vector<matrix> gamma_ij;
		double logprob, oldlogprob, thresh;
		int T, N, M, iters, maxiters;

		HMM(matrix in_A, matrix in_B, matrix in_PI, vector<int> in_O_seq); // constructor
		void forward_pass();
		void backward_pass();
		void gamma_pass();
		void Re_estimate();
		void calclogprob();
		void iterate();
};

// FUNCTION DECLARATIONS AND EXPLANATIONS  HERE

// prints matrix to terminal
void dispmat(matrix inmat);

// Returns transpose of given matrix
matrix transpose (matrix inmat);

// Converts string from input file into a matrix
matrix str2mat (string line);

// Converts string from input file into vector of ints containing the emmission sequence
vector<int> str2seq (string line);

// Converts matrix into string for output format (opposite of str2mat)
string mat2str (matrix inmat);

// creates a row stochastic matrix of NxM with uniform if uni = 1, otherwise random
matrix initmat (int N, int M, int uni);

// MAIN PROGRAM HERE
int main(void)
{
  	string line;
	matrix _A; //transmission matrix
	matrix _B; //emmission matrix
	matrix _PI; //intial state distribution
	vector<int> _O_seq; //emmission sequence
   	
	getline(cin,line);	
	_A = str2mat(line);
	getline(cin, line);
	_B = str2mat(line);
	getline(cin, line);
	_PI = str2mat(line);
	getline(cin, line);
	_O_seq = str2seq(line);

	HMM model(_A, _B, _PI, _O_seq);
	
	model.iterate();

	cout << mat2str(model.A) << "\n" << mat2str(model.B);
}

// FUNCTIONS BODIES HERE

void dispmat(matrix inmat)
{	
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

// creates a row stochastic matrix of NxM with uniform if uni = 1, otherwise random
matrix initmat (int N, int M, int uni)
{
	if(uni == 1)
	{
		matrix outmat (N, vector<double>(M,1.0/M));
		return outmat;
	}
	matrix outmat (N, vector<double>(M));
	for(int i = 0; i < N; i++)
	{
		double sum = 0;
		for(int j = 0; j < M; j++)
		{
			outmat[i][j] = rand();
			sum+= outmat[i][j];
		}
		for(int j = 0; j < M; j++)
		{
			outmat[i][j] /= sum;
		}
	} 
	return outmat;
}

// class functions here

HMM::HMM(matrix in_A, matrix in_B, matrix in_PI, vector<int> in_O_seq)
{
	maxiters = 50000;
	iters = 0;
	thresh = pow(2,-20);
	oldlogprob = -DBL_MAX;
	logprob = -DBL_MAX;
	T = in_O_seq.size();
	N = in_A.size();
	M = in_B[0].size();
	A = in_A;
	B = in_B;
	PI = in_PI;
	O_seq = in_O_seq;
	vector<double> storeprob;
}

void HMM::forward_pass()
{
	matrix mt_alph (T, vector<double>(N)); 
	vector<double> mt_c(T);
	alpha = mt_alph;
	c = mt_c;

	// compute alpha0(i)
	for (int i = 0; i < N; i++)	// from 0 to N-1
	{
		alpha[0][i] = PI[0][i]*B[i][O_seq[0]];
		c[0] += alpha[0][i];
	}

	// scale alpha0(i)
	c[0] = 1.0/c[0];
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
			alpha[t][i] *= B[i][O_seq[t]];
			c[t] += alpha[t][i];
		}
		// scale alphat(i)
		c[t] = 1.0/c[t];
		for (int i = 0; i < N; i++)
		{
			alpha[t][i] *= c[t];
		}
	}
}

void HMM::backward_pass()
{
	// creates beta matrix size TxN filled with zeros
	matrix mt_bet (T, vector<double>(N)); 
	beta = mt_bet;

	// beta_t-1(i) = 1 scaled by C_t-1
	for (int i = 0; i < N; i ++) beta[T-1][i] = c[T-1];

	// backward pass
	for (int t = T-2; t >= 0; t--) // from T-2 to 0
	{
		for (int i = 0; i < N; i ++) // from 0 to N-1
		{
			for (int j = 0; j < N; j++) // from 0 to N-1
			{
				beta[t][i] += A[i][j]*B[j][O_seq[t+1]]*beta[t+1][j];
			}
			beta[t][i] *= c[t];
		}
	}
}

void HMM::gamma_pass()
{
	// creates necessary matrices filled with zeros
	matrix mt_gi (T, vector<double>(N)); 
	vector<matrix> mt_gij (T-1, matrix(N, vector<double>(N)));
	gamma_i = mt_gi;
	gamma_ij = mt_gij;

	// calculates gamma_ij amd gamma_i for t=0 to T-2, making them both length T-1
	for (int t = 0; t < T-1; t++) // from 0 to T-1
	{
		// // Fromw hat I can tell this denom is actually always 1, so this does nothing.
		// double denom = 0;
		// for (int i = 0; i < N; i++) // from 0 to N-1
		// {
		// 	for(int j = 0; j < N; j++)
		// 	{
		// 		denom += alpha[t][i]*A[i][j]*B[j][O_seq[t+1]]*beta[t+1][j];
		// 	}
		// }
		for (int i = 0; i < N; i++) // from 0 to N-1
		{
			for (int j = 0; j < N; j++)
			{
				gamma_ij[t][i][j] = (alpha[t][i]*A[i][j]*B[j][O_seq[t+1]]*beta[t+1][j]);//denom;
				gamma_i[t][i] += gamma_ij[t][i][j]; 
			}
		}
	}
	//special case for gammaT-1(i)
	double denom = 0;
	for (int i = 0; i < N; i++)
	{
		denom += alpha[T-1][i];
	}
	for (int i = 0; i < N; i++)
	{
		gamma_i[T-1][i] = alpha[T-1][i]/denom;
	}
}

void HMM::Re_estimate()
{	
	matrix mt_A (N, vector<double>(N));
	matrix mt_B (N, vector<double>(M));
	matrix mt_PI (1, vector<double>(N));
	A = mt_A;
	B = mt_B;
	PI = mt_PI;

	//re-estimate PI
	for(int i = 0; i < N; i++) PI[0][i] = gamma_i[0][i]; // from O to N-1

	//re-estimate A
	for(int i = 0; i < N; i++)	// from 0 to N-1
	{
		for(int j = 0; j < N; j++) // from 0 to N-1
		{
			double numer = 0, denom = 0;
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
		for(int j = 0; j < M; j++) // from 0 to M-1
		{
			double numer = 0, denom = 0;
			for(int t = 0; t < T; t++) // from 0 to T-1
			{
				if(O_seq[t] == j) {numer+= gamma_i[t][i];}
				denom += gamma_i[t][i];
			}
			B[i][j] = numer/denom;
		}
	}
}

void HMM::calclogprob()
{
	logprob = 0;
	for(int i = 0; i < T; i++) logprob+= log(c[i]);
	logprob = -logprob;
}

void HMM::iterate()
{    
	while((iters < maxiters) && (logprob - thresh >= oldlogprob))
	{	
		oldlogprob = logprob;
		forward_pass();
		backward_pass();
		gamma_pass();
		Re_estimate();
		calclogprob();
		iters ++;
		storeprob.push_back(logprob);
	}
}

