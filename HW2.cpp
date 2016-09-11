/*Compiler settings

For C++, Kattis uses gcc version g++ (Ubuntu 5.4.0-6ubuntu1~16.04.2) 5.4.0 20160609 with the following flags: -g -O2 -static -std=gnu++11 {files}.
System libraries

You are allowed to use all standard libraries included with C++.*/

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

typedef vector< vector<double> > matrix; // matrices are all vectors of vectors of doubles

// FUNCTION DECLARATIONS AND EXPLANATIONS  HERE

// Displays the contents of a matrix to the terminal, columns seperated by tab, rows seperated by newline
void dispmat(matrix inmat); 

// Multiplies two matrices. If there is a dimension mismatch, error is printed to terminal
matrix matmul (matrix mat_a, matrix mat_b);

// Performs dot multiplication of 1 row and 1 column matrix
matrix vecdot (matrix rowmat, matrix colmat);

// returns transpose of given matrix
matrix transpose (matrix inmat);

// Converts string from input file into a matrix
matrix str2mat (string line);

// converts string from input file into vector of ints containing the emmission sequence
vector<int> str2seq (string line);

// Converts matrix into string for output format (opposite of str2mat)
matrix str2mat (string line);

// Returns column matrix which is the given column of the input matrix
matrix matcol (matrix inmat, int col);

// Returns probability of getting first observation of the emmission sequence
matrix alpha1 (matrix B, matrix PI, vector<int> O_seq);

// Does 1 time step for the alpha
matrix alphan (matrix A, matrix B, matrix prev_alpha, int O_t);

// Sums each entry of the alpha matrix
double alphasum (matrix alpha);

// Returns vector of alpha values for each observation and timestep
matrix alpha_forward(matrix A, matrix B, matrix PI, vector<int> O_seq);


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

	matrix outmat = alpha_forward(A, B, PI, O_seq);
	cout << outmat[0][O_seq.size()-1];

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

matrix vecdot (matrix mat_row, matrix mat_col)
{	
	if(mat_row[0].size() != mat_col.size())
	{
		cout << "dimension mismatch\n";
	}
	else
	{
		matrix outmat (1, vector<double> (mat_col.size()) );
		for(int i = 0; i < mat_col.size(); i++)
		{
			outmat[0][i] = mat_row[0][i]*mat_col[i][0];
		}
		return outmat;
	}
}

matrix transpose (matrix inmat)
{
	matrix outmat (inmat[0].size(), vector<double>(inmat.size()));
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
	matrix outmat (1, vector<double> (PI[0].size()) );
	outmat = vecdot(PI,matcol(B,O_seq[0]));
	return outmat;
}

matrix alphan (matrix A, matrix B, matrix prev_alpha, int O_t)
{
	matrix new_alpha (1, vector<double> (prev_alpha[0].size()) );
	matrix alpha_A = matmul(prev_alpha, A);
	matrix B_Ot = matcol(B, O_t);
	new_alpha = vecdot(alpha_A, B_Ot);
	return new_alpha;
}

double alphasum (matrix alpha)
{
	double sum = 0;
	for (int i = 0; i < alpha[0].size(); i++)
	{
		sum += alpha[0][i];
	}
	return sum;
}

matrix alpha_forward(matrix A, matrix B, matrix PI, vector<int> O_seq)
{
	matrix outmat (1, vector<double> (O_seq.size()) );
	matrix prev_alpha = alpha1(B, PI, O_seq);
	outmat[0][0] = alphasum(prev_alpha);

	for (int i = 1; i < O_seq.size(); i++)
	{
		prev_alpha = alphan(A, B, prev_alpha, O_seq[i]);
		outmat[0][i] = alphasum(prev_alpha);
	}
	return outmat;
}