/*Compiler settings

For C++, Kattis uses gcc version g++ (Ubuntu 5.4.0-6ubuntu1~16.04.2) 5.4.0 20160609 with the following flags: -g -O2 -static -std=gnu++11 {files}.
System libraries

You are allowed to use all standard libraries included with C++.*/

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>

using namespace std;

typedef vector< vector<double> > matrix; // matrices are all vectors of vectors of doubles

struct _2mat
{
	matrix val;
	matrix idx;
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
matrix str2mat (string line);

// Returns column matrix which is the given column of the input matrix
matrix matcol (matrix inmat, int col);

// Returns probability of getting first observation of the emmission sequence
matrix alpha1 (matrix B, matrix PI, vector<int> O_seq);

// Does 1 time step for the alpha, returns matrix of alpha values (one for each state of that timestep). Size 1 x size_of_A
matrix alphan (matrix A, matrix B, matrix prev_alpha, int O_t);

// Sums entries of a row matrix size (1xN)
double rowsum (matrix alpha);

// Returns struct containing a matrix (row vector) of alpha values (probability of given observation occuring at that timestep) and a matrix (row vector) of their maximum index (most likely state at that timestep)
vector<int> alpha_forward(matrix A, matrix B, matrix PI, vector<int> O_seq);

// Prints state sequence to cout as required by Kattis
void printseq(vector<int> state_seq);

// Returns normalized row vector given input row vector and normalisation factor (use rowsum function)
matrix normalize(matrix inmat, double normfactor);

// Returns max index of row vector matrix size (1xN)
int maxidxmat(matrix inmat);

// Returns max value of row vector matrix size (1xN)
double maxvalmat(matrix inmat);

// Populates matrix of size rowxcol with double val;
matrix repval(double val, int row, int col);

// Does 1 pass of viterbi and returns struct containing a matrix with the delta values for that timestep, and another for their indexes for that timestep
_2mat viterbin(matrix A, matrix B, matrix prev_delta, int O_t);

// Calculates final delta value (TxN) and delta index matrix (T-1 x N)for all timesteps
_2mat viterbimat(matrix A, matrix B, matrix PI, vector<int> O_seq);

// Given the final delta value and index matrices, calculates the sequence
vector<int> calviterbiseq(matrix val, matrix idx);

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

	_2mat out;
	out = viterbimat(A, B, PI, O_seq);
	vector<int> revans = calviterbiseq(out.val, out.idx);

	for (int i = 1; i < revans.size()+1; i++)
	{
		cout << revans[revans.size()-i] << " ";
	}

	return 0; 
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
	matrix outmat (1, vector<double> (PI[0].size()) ); // creates empty matrix for alpha at first timestep
	outmat = vecdot(PI,matcol(B,O_seq[0])); // multiplies each element of PI by each element of the column of B corresponding to the first emmission. Returns as row vector.
	return outmat;
}

matrix alphan (matrix A, matrix B, matrix prev_alpha, int O_t)
{
	matrix new_alpha (1, vector<double> (prev_alpha[0].size()) ); // creates empty matrix for alpha
	matrix alpha_A = matmul(prev_alpha, A); 
	matrix B_Ot = matcol(B, O_t); // returns column of B corresponding to emmission given.
	new_alpha = vecdot(alpha_A, B_Ot); 
	return new_alpha;
}

double rowsum (matrix alpha)
{
	double sum = 0;
	for (int i = 0; i < alpha[0].size(); i++)
	{
		sum += alpha[0][i];
	}
	return sum;
}

vector<int> alpha_forward (matrix A, matrix B, matrix PI, vector<int> O_seq)
{
	vector<int> maxIDXvec;	// creates vector to hold indices containing state sequence

	matrix prev_alpha = alpha1(B, PI, O_seq); // calculates alpha at first timestep
	double normfac = rowsum(prev_alpha);	// sums values of alpha at first timestep
	prev_alpha = normalize(prev_alpha, normfac);	// normalises alpha at first timestep
	maxIDXvec.push_back(maxidxmat(prev_alpha));		// adds index of highest value of alpha at first timstep to state sequence vector
	
	// Does the same as above for timesteps 2:T
	for (int i = 1; i < O_seq.size(); i++)
	{
		prev_alpha = alphan(A, B, prev_alpha, O_seq[i]);
		normfac = rowsum(prev_alpha);
		prev_alpha = normalize(prev_alpha, normfac);
		maxIDXvec.push_back(maxidxmat(prev_alpha));
	}
	return maxIDXvec;
}

void printseq(vector<int> state_seq)
{
	for(int i = 0; i < state_seq.size(); i++) cout << state_seq[i] << " ";
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

int maxidxmat(matrix inmat)
{
	return(distance(inmat[0].begin(), max_element(inmat[0].begin(), inmat[0].end())));
}

double maxvalmat(matrix inmat)
{
	return *max_element(inmat[0].begin(), inmat[0].end());
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

_2mat viterbin(matrix A, matrix B, matrix prev_delta, int O_t)
{
	_2mat out;
	matrix _val (1, vector<double> (prev_delta[0].size()) );
	matrix _idx (1, vector<double> (prev_delta[0].size()) );
	matrix temp (1, vector<double> (prev_delta[0].size()) );
	matrix bcol;
	for (int i = 0; i < prev_delta[0].size(); i++)
	{
		bcol = repval(B[i][O_t],prev_delta[0].size(),1); // creates column vector populated with the value bi(Ot)
		temp = vecdot(prev_delta,matcol(A,i));	
		temp = vecdot(temp, bcol);
		_val[0][i] = maxvalmat(temp);
		_idx[0][i] = maxidxmat(temp); 
	}
	out.val = _val;
	out.idx = _idx;
	return out;
}

_2mat viterbimat(matrix A, matrix B, matrix PI, vector<int> O_seq)
{	
	matrix prev_delta = alpha1 (B, PI, O_seq); // calculates first delta
	_2mat out;	// creates matrix for deltas and delta indexes
	_2mat final; // creates matrix for final full delta and delta index matrices
	matrix idxsize (O_seq.size()-1, vector<double> (prev_delta.size()) ); // creates empty delta idx matrix of correct size
	matrix valsize (O_seq.size(), vector<double> (prev_delta.size()) ); // creates empty delta value matrix of correct size
	final.val = valsize; //makes the final delta value matrix above size (TxN)
	final.idx = idxsize; //makes the final delta index matrix (TxN)

	final.val[0] = prev_delta[0]; // adds first delta values to delta matrix

	for (int i = 1; i < O_seq.size(); i++) // for each timestep
	{
		out = viterbin(A, B, prev_delta, O_seq[i]); // calulates delta value and delta index matrix (size 1xN) for each timestep
		prev_delta = out.val;
		final.val[i] = prev_delta[0];
		final.idx[i-1] = out.idx[0];
	}

	return final;
}

vector<int> calviterbiseq(matrix val, matrix idx)
{
	vector<int> out;
	int sz = val[0].size();
	matrix temp (1, vector<double> (sz) );

	for (int i = 0; i < sz; i++)
	{
		temp[0][i] = val[sz-1][i];
	}
	int Tfinal = maxidxmat(temp);
	out.push_back(Tfinal);
	int next_state = Tfinal;

	for(int i = 1; i < idx.size() + 1; i++)
	{
		next_state = idx[idx.size()-i][next_state];
		out.push_back(next_state);
	}
	return out;
}