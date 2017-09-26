// HW10.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <vector>
#include "funcs.h"
#include "solveADAMSBASHFORTH.h"
#include "solveADI.h"
#include "solveFuncs.h"
#include "tridiagonal_solver.h"
#include "timestep.h"
#include "funcs.h"

using namespace std;

void recordToFile(ofstream &File, vector< vector<double> > &f, int M, int N)
{
	for (int j = 0; j<N; j++)
	{
		for (int i = 0; i<M; i++)
		{
			File << fixed << setprecision(20) << f[i][j]; // Print to .txt
			if (i != M - 1) { File << "\t"; }
		}
		File << endl; // New Line
	}
}

int main()
{

	ofstream File1;	File1.open("u_data.txt"); // File for u
	ofstream File2;	File2.open("v_data.txt"); // File for v
	ofstream File3;	File3.open("Y_data.txt"); // File for Y

	ofstream File4;	File4.open("uM_data.txt"); // File for u
	ofstream File5;	File5.open("vM_data.txt"); // File for v
	ofstream File6;	File6.open("YM_data.txt"); // File for Y


	const int M = 100;
	const int N = 20;

	const double CFL = 0.85;
	const double alpha = 0.00001;
	const double visc = 200;

	vector< vector<double> > tSlits = { { 1.3, 0.35, 1 },
										{ 3.8, 0.15, 0 } };
	vector< vector<double> > bSlits = { { 0.6, 0.35, 0 },
										{ 1.8, 0.15, 1 } };
	int numSlits = 2;

	int gc = 1;  // Number of ghost cells

	double xLength = 5;
	double yLength = 1;

	double dx = xLength / M; // X grid spacing
	double dy = yLength / N; // Y grid spacing

	double x0 = -dx / 2;
	double x0S = 0;
	double x[M + 2]; // M+2 for cell centered
	double xS[M + 1]; // M+1 for staggered
	double y0 = -dy / 2;
	double y0S = 0;
	double y[N + 2]; // N+2 for cell centered
	double yS[N + 1]; // N+1 for staggered

	vector< vector<double> > u(M + 1, vector<double>(N + 2, 0));
	vector< vector<double> > v(M + 2, vector<double>(N + 1, 0));
	vector< vector<double> > u0(M + 1, vector<double>(N + 2, 0));
	vector< vector<double> > v0(M + 2, vector<double>(N + 1, 0));
	vector< vector<double> > Y(M + 6, vector<double>(N + 6, 0));


	for (int i = 0; i<M + 2; i++) // x centered
	{
		x[i] = x0 + dx*i;
	}

	for (int i = 0; i<N + 2; i++) // y centered
	{
		y[i] = y0 + dy*i;
	}

	for (int i = 0; i<M + 1; i++) // x staggered
	{
		xS[i] = x0S + dx*i;
	}

	for (int i = 0; i<N + 1; i++) // y staggered
	{
		yS[i] = y0S + dy*i;
	}

	uBoundary(u, xS, y, M + 1, N + 2);
	vBoundary(v, x, yS, M + 2, N + 1, tSlits, bSlits, numSlits);
	YBoundary(Y, x, y, M + 2*gc, N + 2 * gc, tSlits, bSlits, numSlits);


	double tmin = 0; // Start of time
	double tmax = 3;
	double dt = timestep(u, v, CFL, dx, dy, M, N);

	double recordTime[5] = { 0.01,0.02,0.03,0.05 }; // Times to write to file
	bool recordZero = true; // Record initial


	bool recordStep = false;
	double t = tmin;
	int record = 0; // Start at 0
	double dtOrig = dt;
	double fps = 10; // Frames for movie
	int fcount = 1; // For movie

	if (recordZero)
	{
		recordToFile(File1, u, M + 1, N + 2); // Record u
		recordToFile(File2, v, M + 2, N + 1); // Record v
		recordToFile(File3, Y, M + 2, N + 2); // Record Y
	}

	// Loop through time
	while (t < tmax)
	{

		dt = timestep(u, v, CFL, dx, dy, M, N);

		/*if (t + dt >= recordTime[record])  // Check for record
		{
		dtOrig = dt;
		dt = recordTime[record] - t;
		if (dt < 0){ break; }
		record++;
		recordStep = true;
		}*/

		solveY2D(Y, u, v, x, y, xS, yS, dx, dy, dt, t, M, N, tSlits, bSlits, numSlits);
		YBoundary(Y, x, y, M + 2 * gc, N + 2 * gc, tSlits, bSlits, numSlits);


		solveBurgers2D(u, v, u0, v0, x, y, dx, dy, dt, t, M, N, tSlits, bSlits, numSlits);
		uBoundary(u, xS, y, M + 1, N + 2);
		vBoundary(v, x, yS, M + 2, N + 1, tSlits, bSlits, numSlits);


		if (recordStep) // Recorded timestep
		{
			recordToFile(File1, u, M + 1, N + 2); // Record u
			recordToFile(File2, v, M + 2, N + 1); // Record v
			recordToFile(File3, Y, M + 2, N + 2); // Record Y

			dt = dtOrig; // Reset to normal step
			recordStep = false;
		}

		if (t - fcount / fps > 0)
		{
			recordToFile(File4, u, M + 1, N + 2); // Record u
			recordToFile(File5, v, M + 2, N + 1); // Record v
			recordToFile(File6, Y, M + 2, N + 2); // Record Y
			fcount++;
		}

		t += dt;
		cout << t << endl;

	}

	File1.close();
	File2.close();
	File3.close();

	system("pause");

}