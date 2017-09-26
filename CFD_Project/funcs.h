#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#pragma once

using namespace std;

double max(vector<double> array, int arrlen) // Function to compute maximum
{
	double prospect = array[0];
	for (int i = 1; i<arrlen; i++)
	{
		if (array[i] > prospect) prospect = array[i];
	}
	return prospect;
}

double min(vector<double> array, int arrlen) // Function to compute maximum
{
	double prospect = array[0];
	for (int i = 1; i<arrlen; i++)
	{
		if (array[i] < prospect) prospect = array[i];
	}
	return prospect;
}

double scalarMin(double a, double b) // Function to compute maximum
{
	if (a > b) return a;
	else return b;
}

double max2d(vector< vector<double> > &array, int arrlen1, int arrlen2) // Function to compute maximum
{
	double prospect = array[0][0];
	for (int i = 0; i<arrlen1; i++)
	{
		for (int j = 0; j<arrlen2; j++)
		{
			if (array[i][j] > prospect) prospect = array[i][j];
		}
	}
	return prospect;
}

double min2d(vector< vector<double> > &array, int arrlen1, int arrlen2) // Function to compute maximum
{
	double prospect = array[0][0];
	for (int i = 0; i<arrlen1; i++)
	{
		for (int j = 0; j<arrlen2; j++)
		{
			if (array[i][j] < prospect) prospect = array[i][j];
		}
	}
	return prospect;
}

void multArrays(vector< vector<double> > &array, vector< vector<double> > &newArray,
	double mult, int M, int N) // Function to multiple array
{
	for (int i = 0; i<M; i++)
	{
		for (int j = 0; j<N; j++)
		{
			newArray[i][j] = array[i][j] * mult;
		}
	}
}

void addArrays(vector< vector<double> > &array1, vector< vector<double> > &array2, vector< vector<double> > &newArray,
	int M, int N) // Function to add 2 arrays
{
	for (int i = 0; i<M; i++)
	{
		for (int j = 0; j<N; j++)
		{
			newArray[i][j] = array1[i][j] + array2[i][j];
		}
	}
}

void setEqual(vector< vector<double> > &array1, vector< vector<double> > &array2,
	int M, int N) // Function to add 2 arrays
{
	for (int i = 0; i<M; i++)
	{
		for (int j = 0; j<N; j++)
		{
			array2[i][j] = array1[i][j];
		}
	}
}


void centerU(vector< vector<double> > &u, vector< vector<double> > &uCC, int M, int N, int gc) // Function to add 2 arrays
{
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			uCC[i+gc][j+gc] = (u[i+1][j+1] + u[i][j+1]) / 2.0;
		}
	}
}

void centerV(vector< vector<double> > &v, vector< vector<double> > &vCC, int M, int N, int gc) // Function to add 2 arrays
{
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			vCC[i+gc][j+gc] = (v[i+1][j+1] + v[i+1][j]) / 2.0;
		}
	}
}

double maxabs2d(vector< vector<double> > &array,
	int M, int N)
{
	double prospect = array[0][0];

	for (int i = 0; i<M; i++)
	{
		for (int j = 0; j<N; j++)
		{
			if (array[i][j] < 0)
			{
				if (array[i][j] * -1 > prospect) prospect = array[i][j] * -1.0;
			}
			else
			{
				if (array[i][j] > prospect) prospect = array[i][j];
			}
		}
	}
	return prospect;
}
