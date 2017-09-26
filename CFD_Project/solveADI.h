#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <vector>
#include "tridiagonal_solver.h"
#include "solveBoundary.h"
#pragma once

using namespace std;

void solveADI(vector< vector<double> > &inp, vector< vector<double> > &source,
	int M, int N, double x[], double y[], double dx, double dy, double dt, int gc, int dir,
	vector< vector<double> > &tSlits, vector< vector<double> > &bSlits, int numSlits)
{

	double d1 = dx / 2;
	double d2 = dy / 2;

	vector< vector<double> > a_arr(M, vector<double>(N, 0));
	vector< vector<double> > b_arr(M, vector<double>(N, 0));
	vector< vector<double> > c_arr(M, vector<double>(N, 0));
	vector< vector<double> > d_arr(M, vector<double>(N, 0));

	vector<double> a_1d(M*N, 0);
	vector<double> b_1d(M*N, 0);
	vector<double> c_1d(M*N, 0);
	vector<double> d_1d(M*N, 0);

	vector< vector<double> > f(M, vector<double>(N, 0)); // Temporary array

	for (int i = 0; i<M; i++) // Borrow input
	{
		for (int j = 0; j<N; j++)
		{
			f[i][j] = inp[i + gc][j + gc];
		}
	}

	double BC;

	// STEP 1
	{
		// Interior
		for (int i = 1; i<M - 1; i++)
		{
			for (int j = 1; j<N - 1; j++)
			{
				a_arr[i][j] = -d1;
				b_arr[i][j] = 1 + 2 * d1;
				c_arr[i][j] = -d1;
				d_arr[i][j] = d2*f[i][j + 1] + (1 - 2 * d2)*f[i][j] + d2*f[i][j - 1] + dt*source[i][j] / 2;
			}
		}

		// Boundaries
		for (int i = 0; i<M; i++)
		{
			a_arr[i][0] = -d1; // Bottom
			b_arr[i][0] = 1 + 2 * d1; // Bottom
			c_arr[i][0] = -d1; // Bottom
			d_arr[i][0] = d2*f[i][1] + (1 - 2 * d2)*f[i][0] + d2*inp[i][0] + dt*source[i][0] / 2; // Bottom
			a_arr[i][N - 1] = -d1; // Top
			b_arr[i][N - 1] = 1 + 2 * d1; // Top
			c_arr[i][N - 1] = -d1; // Top
			d_arr[i][N - 1] = d2*inp[i][N] + (1 - 2 * d2)*f[i][N - 1] + d2*f[i][N - 2] + dt*source[i][N - 1] / 2; // Top
		}

		for (int j = 1; j<N - 1; j++)
		{
			a_arr[0][j] = 0; // Left
			b_arr[0][j] = 1 + 2 * d1; // Left
			c_arr[0][j] = -d1; // Left
			d_arr[0][j] = d2*f[0][j + 1] + (1 - 2 * d2)*f[0][j] + d2*f[0][j - 1] + d1*inp[0][j] + dt*source[0][j] / 2; // Left
			a_arr[M - 1][j] = -d1; // Right
			b_arr[M - 1][j] = 1 + 1 * d1; // Right
			c_arr[M - 1][j] = 0; // Right
			d_arr[M - 1][j] = d2*f[M - 1][j + 1] + (1 - 2 * d2)*f[M - 1][j] + d2*f[M - 1][j - 1] + d1*inp[M][j] + dt*source[M - 1][j] / 2; // Right
		}


		// Load
		int ind = 0;
		for (int j = 0; j<N; j++)
		{
			for (int i = 0; i<M; i++)
			{
				a_1d[ind] = a_arr[i][j];
				b_1d[ind] = b_arr[i][j];
				c_1d[ind] = c_arr[i][j];
				d_1d[ind++] = d_arr[i][j];
			}
		}

		tridiagonal_solver(a_1d, b_1d, c_1d, d_1d, M*N); // Solve tri-diagonal

														 // Unload
		ind = 0;
		for (int j = 0; j<N; j++)
		{
			for (int i = 0; i<M; i++)
			{
				f[i][j] = d_1d[ind++];
			}
		}

		if (dir == 0) uBoundary(f, x, y, M, N);
		else if (dir == 1) vBoundary(f, x, y, M, N, tSlits, bSlits, numSlits);
		else YBoundary(f, x, y, M, N, tSlits, bSlits, numSlits);

	}


	// STEP 2
	{
		// Interior
		for (int i = 1; i<M - 1; i++)
		{
			for (int j = 1; j<N - 1; j++)
			{
				a_arr[i][j] = -d2;
				b_arr[i][j] = 1 + 2 * d2;
				c_arr[i][j] = -d2;
				d_arr[i][j] = d1*f[i + 1][j] + (1 - 2 * d1)*f[i][j] + d1*f[i - 1][j] + dt*source[i][j] / 2;
			}
		}

		// Boundaries
		for (int j = 0; j<N; j++)
		{
			a_arr[0][j] = -d2; // Left
			b_arr[0][j] = 1 + 2 * d2; // Left
			c_arr[0][j] = -d2; // Left
			d_arr[0][j] = d1*f[1][j] + (1 - 1 * d1)*f[0][j] + d2*inp[0][j] + dt*source[0][j] / 2; // Left
			a_arr[M - 1][j] = -d2; // Right
			b_arr[M - 1][j] = 1 + 2 * d2; // Right
			c_arr[M - 1][j] = -d2; // Right
			d_arr[M - 1][j] = d2*inp[M][j] + d1*f[M - 2][j] + (1 - 1 * d1)*f[M - 1][j] + dt*source[M - 1][j] / 2; // Right
		}

		for (int i = 1; i<M - 1; i++)
		{
			a_arr[i][0] = 0; // Bottom
			b_arr[i][0] = 1 + 1 * d2; // Bottom
			c_arr[i][0] = -d2; // Bottom
			d_arr[i][0] = d1*f[i - 1][0] + (1 - 2 * d1)*f[i][0] + d1*f[i + 1][0] + d2*inp[i][0] + dt*source[i][0] / 2; // Bottom
			a_arr[i][N - 1] = -d2; // Top
			b_arr[i][N - 1] = 1 + 1 * d2; // Top
			c_arr[i][N - 1] = 0; // Top
			d_arr[i][N - 1] = d1*f[i - 1][N - 1] + (1 - 2 * d1)*f[i][N - 1] + d1*f[i + 1][N - 1] + d2*inp[i][N] + dt*source[i][N - 1] / 2; // Top
		}

		// Load
		int ind = 0;
		for (int i = 0; i<M; i++)
		{
			for (int j = 0; j<N; j++)
			{
				a_1d[ind] = a_arr[i][j];
				b_1d[ind] = b_arr[i][j];
				c_1d[ind] = c_arr[i][j];
				d_1d[ind++] = d_arr[i][j];
			}
		}

		tridiagonal_solver(a_1d, b_1d, c_1d, d_1d, M*N); // Solve tri-diagonal

														 // Unload
		ind = 0;
		for (int i = 0; i<M; i++)
		{
			for (int j = 0; j<N; j++)
			{
				f[i][j] = d_1d[ind++];
			}
		}

		if (dir == 0) uBoundary(f, x, y, M, N);
		else if (dir == 1) vBoundary(f, x, y, M, N, tSlits, bSlits, numSlits);
		else YBoundary(f, x, y, M, N, tSlits, bSlits, numSlits);

	}

	// Set output	
	for (int i = 0; i<M; i++)
	{
		for (int j = 0; j<N; j++)
		{
			inp[i + gc][j + gc] = f[i][j];
		}
	}

}
