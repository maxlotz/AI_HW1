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

using namespace std;

typedef vector< vector<double> > matrix; // matrices are all vectors of vectors of doubles

struct s1
{
	matrix mat;
	vector<double> vec;
};

struct s2
{
	vector<matrix> mat3d; // essentially a 3d array, its a vector<vector<vector<double>>>
	matrix mat;
};

// FUNCTION DECLARATIONS AND EXPLANATIONS  HERE

// Displays the contents of a matrix to the terminal, columns seperated by tab, rows seperated by newline
void dispmat(matrix inmat); 

// Multiplies two matrices. If there is a dimension mismatch, error is printed to terminal
matrix matmul (matrix mat_a, matrix mat_b);

// Performs dot multiplication of a (NxK) matrix and a 1-column (Kx1) matrix, transposes it and returns a (NxK) matrix
matrix vecdot (matrix inmat, matrix colmat);

// Returns transpose of given matrix
matrix transpose (matrix inmat);

// Converts string from input file into a matrix
matrix str2mat (string line);

// Converts string from input file into vector of ints containing the emmission sequence
vector<int> str2seq (string line);

// Converts matrix into string for output format (opposite of str2mat)
string mat2str (matrix inmat);

// Returns column matrix which is the given column of the input matrix
matrix matcol (matrix inmat, int col);

// Returns probability of getting first observation of the emmission sequence
matrix alpha1 (matrix B, matrix PI, vector<int> O_seq);

// Does 1 time step for the alpha, returns matrix of alpha values (one for each state of that timestep). Size 1 x size_of_A
matrix alphan (matrix A, matrix B, matrix prev_alpha, int O_t);

// Sums entries of a row matrix size (1xN)
double rowsum (matrix alpha);

// Returns struct containing a matrix (row vector) of alpha values (probability of given observation occuring at that timestep) and a matrix (row vector) of their maximum index (most likely state at that timestep)
s1 forward_pass(matrix A, matrix B, matrix PI, vector<int> O_seq);

// Returns normalized row vector given input row vector and normalisation factor (use rowsum function)
matrix normalize(matrix inmat, double normfactor);

// Populates matrix of size rowxcol with double val;
matrix repval(double val, int row, int col);

// Accepts matrix a (TxN) and b (1xN), and populates matrix a row row with values of b, uses row index starting from 0!
matrix popmat(matrix a, matrix b, int row);

// Does 1 timestep for Beta, returns row matrix of Beta values for that timestep
matrix betan(matrix A, matrix B, int O_t, matrix next_beta);

// Performs backward pass and returns beta matrix
matrix backward_pass(vector<double> c, matrix A, matrix B, vector<int> O_seq);

// performs gamma pass, calculates gamma_i (vector) and gamma_ij (matrix) and returns as object of s1
s2 gammas (matrix alphas, matrix betas, matrix A, matrix B, vector<int> O_seq);

// Made timesteps as global variables here.
int T;
int N;
int K;
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

	cout << "\nT is " << T << "\tN is " << N << "\tK is " << K << "\n";

	cout << "\n" << "state transission matrix A size (NxN)" << "\n";
	dispmat(A);		// display A
	cout << "\n";
	cout << "Observation matrix B size (NxK)" << "\n";
	dispmat(B);		// display B
	cout << "\n";
	cout << "Emmission sequence B size (T)" << "\n";
	for (int i = 0; i < T; i++) cout << O_seq[i] << "\t";	// display emmission sequence
	cout << "\n";

	s1 out;
	out = forward_pass(A, B, PI, O_seq);	// calculates alphas
	matrix alphas = out.mat;
	cout << "\n" << "Alpha_t(i) size (TxN)" << "\n";
	dispmat(alphas); // displays alphas
	cout << "\n";

	cout << "Normalisation values c, size (T)" << "\n";
	for(int i = 0; i < T; i++) cout << out.vec[i] << "\t"; // displays normalisation values
	cout << "\n\n";	

	matrix betas = backward_pass(out.vec, A, B, O_seq);		// calculates betas
	cout << "Beta_t(i) size (TxN)" << "\n";
	dispmat(betas);	// displays betas
	cout << "\n";

	s2 out2;

	out2 = gammas(alphas, betas, A, B, O_seq);
	matrix gamma_i = out2.mat;
	vector<matrix> gamma_ij = out2.mat3d;

	cout << "gammat_(i) size (TxN)\n";
	dispmat(gamma_i);
	cout << "\ngammat_(ij) size(T-1xNxN)\n";
	for (int i = 0; i < T-1; i++)
	{
		dispmat(gamma_ij[i]);
		cout << "\n";
	}

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

matrix matmul (matrix mat_a, matrix mat_b)
{
	if (mat_a[0].size() != mat_b.size())
	{
		cout << "dimension mismatch\n";
	}
	else
	{
		matrix mat_c (mat_a.size(),vector<double>(mat_b[0].size()));
		for(int i = 0; i < mat_a.size(); i++)
		{
			for(int j = 0; j < mat_b[0].size(); j++)
			{
				for(int k = 0; k < mat_a[0].size(); k++)
				{
					mat_c[i][j] += mat_a[i][k]*mat_b[k][j]; 
				}
			}

		}   
		return mat_c;
	}
}

matrix vecdot (matrix inmat, matrix mat_col)
{	
	if(inmat[0].size() != mat_col.size())
	{
		cout << "dimension mismatch\n";
	}
	else
	{
		matrix outmat (inmat[0].size(), vector<double> (inmat.size()) );
		for(int i = 0; i < inmat[0].size(); i++)
		{
			for(int j = 0; j < inmat.size(); j++)
			{
				outmat[i][j] = inmat[j][i]*mat_col[i][0];
			}
		}
		return transpose(outmat);
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

matrix matcol (matrix inmat, int col)
{
	matrix outmat (inmat.size(), vector<double> (1));

	for (int i = 0; i < inmat.size(); i++)
	{
		outmat[i][0] = inmat[i][col];
	}
	return outmat;
}

matrix alpha1 (matrix B, matrix PI, vector<int> O_seq)
{
	matrix outmat (1, vector<double>(N)); // creates empty matrix for alpha at first timestep
	outmat = vecdot(PI,matcol(B,O_seq[0])); // multiplies each element of PI by each element of the column of B corresponding to the first emmission. Returns as row vector.
	return outmat;
}

matrix alphan (matrix A, matrix B, matrix prev_alpha, int O_t)
{
	matrix new_alpha (1, vector<double>(N)); // creates empty matrix for alpha
	matrix alpha_A = matmul(prev_alpha, A); 
	matrix B_Ot = matcol(B, O_t); // returns column of B corresponding to emmission given.
	new_alpha = vecdot(alpha_A, B_Ot); 
	return new_alpha;
}

double rowsum (matrix alpha)
{
	double sum = 0;
	for (int i = 0; i < N; i++)
	{
		sum += alpha[0][i];
	}
	return sum;
}

s1 forward_pass (matrix A, matrix B, matrix PI, vector<int> O_seq)
{
	s1 out;
	matrix alphas (T, vector<double>(N)); // creates empty matrix size (TxN) for alphas
	vector<double> c; // creates vector of doubles to hold normalisation factors c

	matrix prev_alpha = alpha1(B, PI, O_seq); // calculates alpha at first timestep
	double normfac = rowsum(prev_alpha);	// sums values of alpha at first timestep
	prev_alpha = normalize(prev_alpha, normfac);	// normalises alpha at first timestep

	alphas = popmat(alphas, prev_alpha, 0);
	c.push_back(normfac);

	// Does the same as above for timesteps 2:T
	for (int i = 1; i < T; i++)
	{
		prev_alpha = alphan(A, B, prev_alpha, O_seq[i]);
		normfac = rowsum(prev_alpha);
		prev_alpha = normalize(prev_alpha, normfac);

		alphas = popmat(alphas, prev_alpha, i);
		c.push_back(normfac);
	}
	out.mat = alphas;
	out.vec = c;
	return out;
}

matrix normalize(matrix inmat, double normfactor)
{
	matrix outmat (1, vector<double> (inmat[0].size()) );
	for(int i = 0; i < inmat[0].size(); i++)
	{
		outmat[0][i] = inmat[0][i]/normfactor;
	}
	return outmat;
}

matrix repval(double val, int row, int col)
{
	matrix outmat (row, vector<double> (col) );
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < col; j++)
		{
			outmat[i][j] = val;
		}
	}
	return outmat;
}

matrix popmat(matrix a, matrix b, int row)
{
	for (int i = 0; i < a[0].size(); i++)
	{
		a[row][i] = b[0][i];
	}
	return a;
}

matrix betan(matrix A, matrix B, int O_t, matrix next_beta)
{
	matrix prev_beta = vecdot(next_beta, matcol(B, O_t));
	prev_beta = matmul(A, transpose(prev_beta));
	return transpose(prev_beta);
}

matrix backward_pass(vector<double> c, matrix A, matrix B, vector<int> O_seq)
{
	matrix betas (T,vector<double>(N));
	matrix beta1 (1, vector<double> (N,1));
	beta1 = normalize(beta1,c[T-1]);
	betas = popmat(betas, beta1, T-1);
	matrix prev_beta = beta1;

	for (int i = T-2; i > -1; i--) //from T-2 to 0
	{
		prev_beta = betan(A, B, O_seq[i+1], prev_beta);
		prev_beta = normalize(prev_beta, c[i]);
		betas = popmat(betas, prev_beta, i);
	}
	return betas;
}

s2 gammas (matrix alphas, matrix betas, matrix A, matrix B, vector<int> O_seq)
{
	s2 out;
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
				gamma_ij[t][i][j] = (alphas[t][i]*A[i][j]*B[j][O_seq[t+1]])/denom;
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
	out.mat = gamma_i;
	return out;
}