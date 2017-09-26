#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include "solveBoundary.h"
#pragma once

using namespace std;

void weno(vector< vector<double> > &Y, vector< vector<double> > &dYdx, vector< vector<double> > &dYdy,
	vector< vector<double> > &u, vector< vector<double> > &v, int M, int N, double dx, double dy, int num_GC)
	
	{
		double a;
		double b;
		double c;
		double d;
		
		double e = 1e-6;
		
		double IS0;
		double IS1;
		double IS2;
		
		double a0;
		double a1;
		double a2;
		
		double w0;
		double w2;
		
		double func;
	
	for (int i=num_GC; i<M+num_GC; i++)	// For x velocity
	{
		for (int j=num_GC; j<N+num_GC; j++)	
		{
			if (u[i][j] >= 0)
			{
				a = ((Y[i-1][j] - Y[i-2][j]) - (Y[i-2][j] - Y[i-3][j]) ) / dx;
				b = ((Y[i][j]   - Y[i-1][j]) - (Y[i-1][j] - Y[i-2][j]) ) / dx;
				c = ((Y[i+1][j] - Y[i][j]  ) - (Y[i][j]   - Y[i-1][j]) ) / dx;
				d = ((Y[i+2][j] - Y[i+1][j]) - (Y[i+1][j] - Y[i][j]  ) ) / dx;
			} else {
				a = ((Y[i+3][j] - Y[i+2][j]) - (Y[i+2][j] - Y[i+1][j]) ) / dx;
				b = ((Y[i+2][j] - Y[i+1][j]) - (Y[i+1][j] - Y[i][j]  ) ) / dx;
				c = ((Y[i+1][j] - Y[i][j]  ) - (Y[i][j]   - Y[i-1][j]) ) / dx;
				d = ((Y[i][j]   - Y[i-1][j]) - (Y[i-1][j] - Y[i-2][j]) ) / dx;
			}
			
			IS0 = 13*(a-b)*(a-b) + 3*(a-3*b)*(a-3*b);
			IS1 = 13*(b-c)*(b-c) + 3*(b+c)*(b+c);
			IS2 = 13*(c-d)*(c-d) + 3*(3*c-d)*(3*c-d);
		
			a0 = 1 / ((e + IS0)*(e + IS0));
			a1 = 6 / ((e + IS1)*(e + IS1));
			a2 = 3 / ((e + IS2)*(e + IS2));
	
			w0 = a0 / (a0+a1+a2);
			w2 = a2 / (a0+a1+a2);
			
			
			func = (1.0/3)*w0*(a-(2*b)+c) + (1.0/6)*(w2-(1.0/2))*(b-(2*c)+d);
			
			if (u[i][j] >= 0)
			{
				dYdx[i][j] =  (1/(12*dx))	*	( -   (Y[i-1][j]-Y[i-2][j])
					     			      + 7.0*(Y[i][j]-Y[i-1][j])
					                      + 7.0*(Y[i+1][j]-Y[i][j])
			                    		  -   (Y[i+2][j]-Y[i+1][j])) - func;
			} else {
				dYdx[i][j] =  (1/(12*dx))	*	( -   (Y[i-1][j]-Y[i-2][j])
					     				  + 7.0*(Y[i][j]-Y[i-1][j])
					     				  + 7.0*(Y[i+1][j]-Y[i][j])
					    				  -   (Y[i+2][j]-Y[i+1][j])) + func;
			}
		}
	}

	for (int i=num_GC; i<M+num_GC; i++)	// For y velocity
	{
		for (int j=num_GC; j<N+num_GC; j++)	
		{
			if (v[i][j] >= 0)
			{
				a = ((Y[i][j-1] - Y[i][j-2]) - (Y[i][j-2] - Y[i][j-3]) ) / dy;
				b = ((Y[i][j]   - Y[i][j-1]) - (Y[i][j-1] - Y[i][j-2]) ) / dy;
				c = ((Y[i][j+1] - Y[i][j]  ) - (Y[i][j]   - Y[i][j-1]) ) / dy;
				d = ((Y[i][j+2] - Y[i][j+1]) - (Y[i][j+1] - Y[i][j]  ) ) / dy;
			} else {
				a = ((Y[i][j+3] - Y[i][j+2]) - (Y[i][j+2] - Y[i][j+1]) ) / dy;
				b = ((Y[i][j+2] - Y[i][j+1]) - (Y[i][j+1] - Y[i][j]  ) ) / dy;
				c = ((Y[i][j+1] - Y[i][j]  ) - (Y[i][j]   - Y[i][j-1]) ) / dy;
				d = ((Y[i][j]   - Y[i][j-1]) - (Y[i][j-1] - Y[i][j-2]) ) / dy;
			}
			
			IS0 = 13*(a-b)*(a-b) + 3*(a-3*b)*(a-3*b);
			IS1 = 13*(b-c)*(b-c) + 3*(b+c)*(b+c);
			IS2 = 13*(c-d)*(c-d) + 3*(3*c-d)*(3*c-d);
		
			a0 = 1 / ((e + IS0)*(e + IS0));
			a1 = 6 / ((e + IS1)*(e + IS1));
			a2 = 3 / ((e + IS2)*(e + IS2));
	
			w0 = a0 / (a0+a1+a2);
			w2 = a2 / (a0+a1+a2);
			
			
			func = (1.0/3)*w0*(a-(2*b)+c) + (1.0/6)*(w2-(1.0/2))*(b-(2.0*c)+d);
			
			if (v[i][j] >= 0)
			{
				dYdy[i][j] =  (1/(12.0*dx))	*	( -   (Y[i][j-1]-Y[i][j-2])
					     			      + 7.0*(Y[i][j]-Y[i][j-1])
					                      + 7.0*(Y[i][j+1]-Y[i][j])
			                    		  -   (Y[i][j+2]-Y[i][j+1])) - func;
			} else {
				dYdy[i][j] =  (1/(12.0*dx))	*	( -   (Y[i][j-1]-Y[i][j-2])
					     				  + 7.0*(Y[i][j]-Y[i][j-1])
					     				  + 7.0*(Y[i][j+1]-Y[i][j])
					    				  -   (Y[i][j+2]-Y[i][j+1])) + func;
			}
		}
	}
}

void solveWENO5RK3( vector< vector<double> > &Y, vector< vector<double> > &u, vector< vector<double> > &v,
					double x[], double y[], double dx, double dy, double dt, double t, int M, int N, int gc,
					vector< vector<double> > &tSlits, vector< vector<double> > &bSlits, int numSlits)
{
	
	vector< vector<double> > dYdx0(M + 2 * gc, vector<double>(M + 2 * gc, 0));
	vector< vector<double> > dYdx1(M + 2 * gc, vector<double>(M + 2 * gc, 0));
	vector< vector<double> > dYdx2(M + 2 * gc, vector<double>(M + 2 * gc, 0));
	
	vector< vector<double> > dYdy0(M + 2 * gc, vector<double>(N + 2 * gc, 0));
	vector< vector<double> > dYdy1(M + 2 * gc, vector<double>(N + 2 * gc, 0));
	vector< vector<double> > dYdy2(M + 2 * gc, vector<double>(N + 2 * gc, 0));
	
	vector< vector<double> > Y1(M + 2 * gc, vector<double>(N + 2 * gc, 0));
	vector< vector<double> > Y2(M + 2 * gc, vector<double>(N + 2 * gc, 0));
	
	double a10 = 1.0;
	double a20 = -3.0/4;
	double a21 = 1.0/4;
	double a30 = -1.0/12;
	double a31 = -1.0/12;
	double a32 = 2.0/3;

	{ // Step 1
	
		weno(Y, dYdx0, dYdy0, u, v, M, N, dx, dy, 3); // Weno
	
		for (int i=0; i<M + gc * 2; i++) // RK
		{
			for (int j=0; j<N + gc * 2; j++) // RK
			{
				
				Y1[i][j] = Y[i][j] - a10*(u[i][j]*dt*dYdx0[i][j] + v[i][j]*dt*dYdy0[i][j]);
			}
		}
		
		YBoundary(Y1, x, y, M + 2 * gc, N + 2 * gc, tSlits, bSlits, numSlits); // Boundary
		
	}
	
	{ // Step 2
		weno(Y1, dYdx1, dYdy1, u, v, M, N, dx, dy, 3); // Weno
		
		for (int i=0; i<M + gc * 2; i++) // RK
		{
			for (int j=0; j<N + gc * 2; j++) // RK
			{
				Y2[i][j] = Y1[i][j] - a20*(u[i][j]*dt*dYdx0[i][j] + v[i][j]*dt*dYdy0[i][j])
								    - a21*(u[i][j]*dt*dYdx1[i][j] + v[i][j]*dt*dYdy1[i][j]);
			}
		}
		
		YBoundary(Y2, x, y, M + 2 * gc, N + 2 * gc, tSlits, bSlits, numSlits); // Boundary
	}
	
	{ // Step 3
		weno(Y2, dYdx2, dYdy2, u, v, M , N, dx, dy, 3); // Weno
		
		for (int i=0; i<M + gc * 2; i++) // RK
		{
			for (int j=0; j<N + gc * 2; j++) // RK
			{
				Y[i][j] = Y2[i][j] - a30*(u[i][j]*dt*dYdx0[i][j] + v[i][j]*dt*dYdy0[i][j])
								   - a31*(u[i][j]*dt*dYdx1[i][j] + v[i][j]*dt*dYdy1[i][j])
				                   - a32*(u[i][j]*dt*dYdx2[i][j] + v[i][j]*dt*dYdy2[i][j]);
			}
		}
		
		YBoundary(Y, x, y, M + 2 * gc, N + 2 * gc, tSlits, bSlits, numSlits); // Boundary
	}
	
}
