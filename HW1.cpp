#include <iostream>
#include <vector>
#include <tuple>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

typedef vector< vector<double> > matrix; // matrices are all vectors of vectors of doubles

// FUNCTION DECLARATIONS AND EXPLANATIONS  HERE

// displays the contents of a matrix to the terminal, columns seperated by tab, rows seperated by newline
void dispmat(vector< vector<double> >); 

//multiplies two matrices, prints to terminal if dimensions aren't consistent
vector< vector<double> > matmul (vector< vector<double> >, vector< vector<double> >);

// creates two matrices for testing
tuple< vector< vector<double> >, vector< vector<double> > > testinit(void);

// Converts string from input file into a matrix
vector< vector<double> > str2mat (string);

// Converts matrix into string for output format (opposite of str2mat)
string mat2str (vector< vector<double> >);

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
		vector< vector<double> > A; //transmission matrix
		vector< vector<double> > B; //emmission matrix
		vector< vector<double> > PI; //intial state distribution
	   	
	   	string line;
	   	string D;
		getline(infile, line);
		A = str2mat(line);
		D = mat2str(A);
		cout << D;

		return 0;
	}
}

// FUNCTIONS BODIES HERE
void dispmat(vector<vector<double>> vec)
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

vector< vector<double> > matmul (vector< vector<double> > matrixa, vector< vector<double> > matrixb)
{
	vector< vector<double> > matrixc (matrixa.size(),vector<double>(matrixb[0].size()));
	if (matrixa[0].size() != matrixb.size())
	{
		cout << "matrix dimensions inconsistent \n";
	}
	else
	{
		for(int i = 0; i < matrixa.size(); i++)
		{
			for(int j = 0; j < matrixb[0].size(); j++)
			{
				for(int k = 0; k < matrixa[0].size(); k++)
				{
					matrixc[i][j] += matrixa[i][k]*matrixb[k][j]; 
				}
			}

		}   
		return matrixc;
	}
}

tuple< vector< vector<double> >, vector< vector<double> > > testinit(void)
{
	vector< vector<double> > matrixa;
	vector<double> a_row1{0.1,0.2};
	vector<double> a_row2{0.3,0.4};
	vector<double> a_row3{0.5,0.6};
	matrixa.push_back(a_row1);
	matrixa.push_back(a_row2);
	matrixa.push_back(a_row3);

	vector< vector<double> > matrixb;
	vector<double> b_row1{0.05, 0.1, 0.15, 0.2, 0.25, 0.3};
	vector<double> b_row2{0.35, 0.4, 0.45, 0.5, 0.55, 0.6};
	matrixb.push_back(b_row1);
	matrixb.push_back(b_row2);

	return make_tuple(matrixa, matrixb);
}

vector< vector<double> > str2mat (string line)
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
	vector< vector<double> > outmat (rows, vector<double> (cols));

	for(int i = 0; i < rows; i++)
	{
		for(int j = 0; j < cols; j++)
		{
			outmat[i][j] = stod(record[i*cols + j + 2]);
		}
	}
	return outmat;
}

string mat2str (vector< vector<double> > inmat)
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
			outline += blank + to_string(inmat[i][j]); 
		}
	}
	return outline;
}