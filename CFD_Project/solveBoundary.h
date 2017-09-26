#pragma once
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>


void uBoundary(vector< vector<double> > &f, double x[], double y[], int M, int N)
{

	// Left
	for (int j = 0; j < N; j++)
	{
		if (y[j] < 0.5) // Bottom inlet
		{
			f[0][j] = -16 * ((y[j] - 0.25)*(y[j] - 0.25)) + 1;
		}
		else // Top inlet
		{
			f[0][j] = -16 * ((y[j] - 0.75)*(y[j] - 0.75)) + 1;
		}
	}

	// Right
	for (int j = 0; j < N; j++)
	{
		f[M - 1][j] = f[M - 2][j];
	}

	// Bottom
	for (int i = 0; i < M; i++)
	{
		f[i][0] = 0 - f[i][1];
	}

	// Top
	for (int i = 0; i < M; i++)
	{
		f[i][N - 1] = 0 - f[i][N - 2];
	}

}

void vBoundary(vector< vector<double> > &f, double x[], double y[], int M, int N, vector< vector<double> > &tSlits, vector< vector<double> > &bSlits, int numSlits)
{

	// Left
	for (int j = 0; j < N; j++)
	{
		f[0][j] = -f[1][j];
	}

	// Right
	for (int j = 0; j < N; j++)
	{
		f[M - 1][j] = -f[M - 2][j];
	}

	// Bottom
	for (int i = 0; i < M; i++)
	{
		for (int s = 0; s < numSlits; s++) // Loop for each slit
		{
			if ((x[i] >= bSlits[s][0] - bSlits[s][1] / 2) && (x[i] <= bSlits[s][0] + bSlits[s][1] / 2))
			{
				f[i][0] = 0.5;
				break;
			}
			else
			{
				f[i][0] = 0;
			}
		}
	}

	// Top
	for (int i = 0; i < M; i++)
	{
		for (int s = 0; s < numSlits; s++) // Loop for each slit
		{
			if ((x[i] >= tSlits[s][0] - tSlits[s][1] / 2) && (x[i] <= tSlits[s][0] + tSlits[s][1] / 2))
			{
				f[i][N - 1] = -0.5;
				break;
			}
			else
			{
				f[i][N - 1] = 0;
			}
		}
	}

}


void YBoundary(vector< vector<double> > &f, double x[], double y[], int M, int N, vector< vector<double> > &tSlits, vector< vector<double> > &bSlits, int numSlits)
{

	// Left
	for (int j = 0; j < N; j++)
	{
		if (y[j] < 0.5) // Bottom inlet
		{
			f[0][j] = 2 - f[5][j];
			f[1][j] = 2 - f[4][j];
			f[2][j] = 2 - f[3][j];
		}
		else // Top inlet
		{
			f[0][j] = 0 - f[5][j];
			f[1][j] = 0 - f[4][j];
			f[2][j] = 0 - f[3][j];
		}
	}

	// Right
	for (int j = 0; j < N; j++)
	{
		f[M - 1][j] = f[M - 6][j];
		f[M - 2][j] = f[M - 5][j];
		f[M - 3][j] = f[M - 4][j];
	}

	// Top
	for (int i = 0; i < M; i++)
	{
		for (int s = 0; s < numSlits; s++) // Loop for each slit
		{
			if ((x[i] >= tSlits[s][0] - tSlits[s][1] / 2) && (x[i] <= tSlits[s][0] + tSlits[s][1] / 2))
			{
				f[i][N - 1] = 2 * tSlits[s][2] - f[i][N - 6];
				f[i][N - 2] = 2 * tSlits[s][2] - f[i][N - 5];
				f[i][N - 3] = 2 * tSlits[s][2] - f[i][N - 4];
				break;
			}
			else
			{
				f[i][N - 1] = f[i][N - 6];
				f[i][N - 2] = f[i][N - 5];
				f[i][N - 3] = f[i][N - 4];
			}
		}
	}

	// Bottom
	for (int i = 0; i < M; i++)
	{
		for (int s = 0; s < numSlits; s++) // Loop for each slit
		{
			if ((x[i] >= bSlits[s][0] - bSlits[s][1] / 2) && (x[i] <= bSlits[s][0] + bSlits[s][1] / 2))
			{
				f[i][0] = 2 * bSlits[s][2] - f[i][5];
				f[i][1] = 2 * bSlits[s][2] - f[i][4];
				f[i][2] = 2 * bSlits[s][2] - f[i][3];
				break;
			}
			else
			{
				f[i][0] = f[i][5];
				f[i][1] = f[i][4];
				f[i][2] = f[i][3];
			}
		}
	}

}

void genBoundary(vector< vector<double> > &f, double x[], double y[], int M, int N, int gc, int dir, vector< vector<double> > &tSlits, vector< vector<double> > &bSlits, int numSlits)
{
	if (dir == 0) uBoundary(f, x, y, M, N);
	else if (dir == 1) vBoundary(f, x, y, M, N, tSlits, bSlits, numSlits);
	else YBoundary(f, x, y, M, N, tSlits, bSlits, numSlits);
}
