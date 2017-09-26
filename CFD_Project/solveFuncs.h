#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include "solveADAMSBASHFORTH.h"
#include "solveWENO5RK3.h"
#include "solveADI.h"
#include "funcs.h"
#include "solveBoundary.h"

using namespace std;

// Function to Solve Crank Nicholson ADI
void solveBurgers2D(	vector< vector<double> > &u,  vector< vector<double> > &v,
 						vector< vector<double> > &u0,  vector< vector<double> > &v0,
						double x[], double y[], double dx, double dy, double dt, double t, int M, int N,
						vector< vector<double> > &tSlits, vector< vector<double> > &bSlits, int numSlits)
{

	int gc = 1;

	vector< vector<double> > source1(M, vector<double>(N+1, 0));
	vector< vector<double> > source2(M+1, vector<double>(N, 0));

	// Get Source
	solveADAMSBASHFORTH(u, v, u0, v0, source1, dx, dy, M, N+1, gc); // x velocity
	solveADAMSBASHFORTH(v, u, u0, v0, source2, dx, dy, M+1, N, gc); // y velocity
	
	// Call ADI
	solveADI(u, source1, M, N+1, x, y, dx, dy, dt, gc, 1, tSlits, bSlits, numSlits); // x velocity
	solveADI(v, source2, M+1, N, x, y, dx, dy, dt, gc, 2, tSlits, bSlits, numSlits); // y velocity
	
	uBoundary(u, x, y, M+1, N+2);
	vBoundary(v, x, y, M+2, N+1, tSlits, bSlits, numSlits);

	setEqual(u, u0, M+1, N+2);
	setEqual(v, v0, M+2, N+1);
	
}


// Function to Solve Y Equation ADI
void solveY2D( 	vector< vector<double> > &Y,  vector< vector<double> > &u,  vector< vector<double> > &v,
				double x[], double y[], double xS[], double yS[], double dx, double dy, double dt, double t, int M, int N,
				vector< vector<double> > &tSlits, vector< vector<double> > &bSlits, int numSlits)
{
	
	int gc = 3;
	
	vector< vector<double> > source(M + 2 * gc, vector<double>(N + 2 * gc, 0));
	
	vector< vector<double> > uCC(M + 2 * gc, vector<double>(N + 2 * gc, 0)); // Cell center the x velocity
	vector< vector<double> > vCC(M + 2 * gc, vector<double>(N + 2 * gc, 0)); // Cell center the y velocity
	
	centerU(u, uCC, M, N, gc);
	centerV(v, vCC, M, N, gc);
	uBoundary(u, xS, y, M + 1, N + 2);
	vBoundary(v, x, yS, M + 2, N + 1, tSlits, bSlits, numSlits);

	for (int i=0; i<M+2; i++)
	{
		for (int j=0; j<N+2; j++)
		{
			source[i][j] = Y[i][j];
		}
	}
	
	// Get Source
	solveWENO5RK3(source, uCC, vCC, x, y, dx, dy, dt, t, M, N, gc, tSlits, bSlits, numSlits); // Get the source

	for (int i=0; i<M; i++)
	{
		for (int j=0; j<N; j++)
		{
			source[i][j] = source[i][j]-Y[i][j];
		}
	}
	
	// Call ADI
	solveADI(Y, source, M-1, N-1, x, y, dx, dy, dt, gc, 3, tSlits, bSlits, numSlits);

	
}
