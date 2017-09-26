#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <vector>
#pragma once

using namespace std;

// Function to get Adams-Bashforth Term
void solveADAMSBASHFORTH( vector< vector<double> > &u,  vector< vector<double> > &v, vector< vector<double> > &uO,  vector< vector<double> > &vO,
					vector< vector<double> > &adams, double dx, double dy, int M, int N, int gc)
{
	for (int i=1; i<M-1; i++)
	{
		for (int j=1; j<N-1; j++)
		{
			double ft = (1/dx) * ( // Calculate first term
								 pow((0.5*( u[i+1][j]-u[i][j] )),2)
								-pow((0.5*( u[i][j]-u[i-1][j] )),2)
								);
	
			double st = (1/dy) * ( // Calculate second term
								 ((0.5*( uO[i][j]-uO[i][j+1] ))*(0.5*( vO[i][j]-vO[i+1][j] )))
								-((0.5*( uO[i][j-1]-uO[i][j] ))*(0.5*( vO[i][j-1]-vO[i+1][j-1] )))
								);
								
			adams[i-1][j-1] = (1.5*ft - 0.5*st); 
		}
	}
}
