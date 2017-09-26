#include <iostream>
#pragma once

using namespace std;


/* Function to solve tri diagonal matrices

a_arr = array of numbers below diagonal
b_arr = array of numbers on diagonal
c_arr = array of numbers above diagonal
a_arr = right side of equation array
arrLength = dimension of matrix

*/

void tridiagonal_solver(vector<double> a_arr, vector<double> b_arr, vector<double> c_arr, vector<double> d_arr, int arrLength)
{

	// Forward Elimination

	for (int i = 1; i<arrLength; i++)
	{
		b_arr[i] = b_arr[i] - c_arr[i - 1] * a_arr[i] / b_arr[i - 1];
		d_arr[i] = d_arr[i] - d_arr[i - 1] * a_arr[i] / b_arr[i - 1];
	}

	// Backward Elimination

	d_arr[arrLength - 1] = d_arr[arrLength - 1] / b_arr[arrLength - 1];

	for (int i = arrLength - 2; i >= 0; i--)
	{
		d_arr[i] = (d_arr[i] - (c_arr[i] * d_arr[i + 1])) / b_arr[i];
	}

}
