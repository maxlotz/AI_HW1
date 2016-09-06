// ALL MATRICES ARE type vector< vector<double> >

#include <iostream>
#include <vector>
#include <tuple>

using namespace std;

// FUNCTION DECLARATIONS AND EXPLANATIONS  HERE

// displays the contents of a matrix to the terminal, columns seperated by tab, rows seperated by newline
void dispmat(vector< vector<double> >); 

//multiplies two matrices, returns -1 if dimensions arent consistent
vector< vector<double> > matmul (vector< vector<double> >, vector< vector<double> >);

// creates two matrices for testing
tuple< vector< vector<double> >, vector< vector<double> > > testinit(void);

// MAIN PROGRAM HERE
int main(void)
{
	vector< vector<double> > matrixa;
	vector< vector<double> > matrixb;
	vector< vector<double> > matrixc;

	tie(matrixa, matrixb) = testinit();
	matrixc = matmul(matrixa, matrixb);

	cout <<"matrix a:\n";
	dispmat(matrixa);
	cout << "\n";
	cout <<"matrix b:\n";
	dispmat(matrixb);
	cout << "\n";
	cout <<"matrix a x matrix b = matrix c:\n";
	dispmat(matrixc);

	return 0;
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