#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

typedef vector< vector<double> > matrix; // matrices are all vectors of vectors of doubles

// FUNCTION DECLARATIONS AND EXPLANATIONS  HERE

// Displays the contents of a matrix to the terminal, columns seperated by tab, rows seperated by newline
void dispmat(matrix); 

// Multiplies two matrices. If matrix dimensions aren't consistent, a matrix populated with 0s will be returned, and an error message will be printed to the terminal.
matrix matmul (matrix, matrix);

// returns transpose of given matrix
matrix transpose (matrix);

// Converts string from input file into a matrix
matrix str2mat (string);

// converts string from input file into vector of ints containing the emmission sequence
vector<int> str2seq (string);

// Converts matrix into string for output format (opposite of str2mat)
string mat2str (matrix);

// MAIN PROGRAM HERE
int main(int argc, char *argv[])
{

    ifstream infile ( argv[1] );
    if ( !infile.is_open() )
    {
    	cout<<"Could not open file\n";
    }
    else 
    {
		matrix A; //transmission matrix
		matrix B; //emmission matrix
		matrix PI; //intial state distribution
	   	
	 	string line;
		getline(infile, line);
		A = str2mat(line);
		getline(infile, line);
		B = str2mat(line);
		getline(infile, line);
		PI = str2mat(line);
		//getline(infile, line);
		//O_seq = str2seq(line); // emmission sequence

		return 0; // answer goes here, eg. return mat2str(answer_matrix);
	}
}

// FUNCTIONS BODIES HERE
void dispmat(matrix vec)
{
	for (int i = 0; i < vec.size(); i++)
	{
	    for (int j = 0; j < vec[i].size(); j++)
	    {
	        cout << vec[i][j] << "\t";
	    }
	    cout << "\n";
	}
}

matrix matmul (matrix mat_a, matrix mat_b)
{
	matrix mat_c (mat_a.size(),vector<double>(mat_b[0].size()));
	if (mat_a[0].size() != mat_b.size())
	{
		cout << "matrix dimensions inconsistent\n";
	}
	else
	{
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

matrix transpose (matrix inmat)
{
	matrix outmat (inmat[0].size(), vector<double>(inmat.size()));
	for (int i = 0; i < inmat.size(); i++)
	{
		for (int j = 0; j < inmat[0].size(); j++)
		{
			outmat[i][j] = inmat[j][i];
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