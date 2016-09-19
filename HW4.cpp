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
#include <float.h>

using namespace std;

typedef vector< vector<double> > matrix; // matrices are all vectors of vectors of doubles

class HMM {
	private:

	public:
		matrix A, B, PI, alpha, beta, gamma_i;
		vector<int> O_seq;
		vector<double> c;
		vector<matrix> gamma_ij;
		double logprob, oldlogprob;
		int T, N, K, iters, maxiters;

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
	cout << "\niteration reached = " << model.iters << "\nlogprob = " << model.logprob << "\noldlogprob = " << model.oldlogprob << "\n\n";
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

// class functions here

HMM::HMM(matrix in_A, matrix in_B, matrix in_PI, vector<int> in_O_seq)
{
	maxiters = 1000;
	iters = 0;
	oldlogprob = -DBL_MAX;
	T = in_O_seq.size();
	N = in_A.size();
	K = in_B[0].size();
	A = in_A;
	B = in_B;
	PI = in_PI;
	O_seq = in_O_seq;
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
}

void HMM::backward_pass()
{
	// creates beta matrix size TxN filled with zeros
	matrix mt_bet (T, vector<double>(N)); 
	beta = mt_bet;

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
		double denom = 0;
		for (int i = 0; i < N; i ++) // from 0 to N-1
		{
			for(int j = 0; j < N; j++)
			{
				denom += alpha[t][i]*A[i][j]*B[j][O_seq[t+1]]*beta[t+1][j];
			}
		}
		for (int i = 0; i < N; i++) // from 0 to N-1
		{
			for (int j = 0; j < N; j++)
			{
				gamma_ij[t][i][j] = (alpha[t][i]*A[i][j]*B[j][O_seq[t+1]]*beta[t+1][j])/denom;
				gamma_i[t][i] += gamma_ij[t][i][j]; 
			}
		}
	}
	// special case for gammaT-1(i)
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
}

void HMM::calclogprob()
{
	logprob = 0;
	for(int i = 0; i < T; i++) logprob+= log(c[i]);
	logprob = -logprob;
}

void HMM::iterate()
{    
	while((iters < maxiters) && (logprob >= oldlogprob))
	{
		iters ++;
		forward_pass();
		backward_pass();
		gamma_pass();
		Re_estimate();
		calclogprob();
		oldlogprob = logprob;
	}
}