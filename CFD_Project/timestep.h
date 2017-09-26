#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <vector>
#include "funcs.h"
#pragma once

double timestep(vector< vector<double> > &u, vector< vector<double> > &v, double CFL, double dx, double dy, int M, int N)
{
	
	double dtu = dx / (2*maxabs2d(u,M+1,N+2)+maxabs2d(v,M+2,N+1)); // Check u
	double dtv = dy / (2*maxabs2d(v,M+1,N+1)+maxabs2d(u,M+1,N+2)); // Check v

	return CFL * scalarMin(dtu,dtv); // Return the minimum
}




