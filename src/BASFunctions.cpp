// I'm just playing around on this to see if I can make the Halton Sequence in C++.
// Next phase is to build this into the Master Sample.
//

#include <iostream>
#include <math.h>
#include <numeric>
#include <vector>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double log_b(double x, int & base)
{
	double tmp;
	tmp = log(x) / log(base);
	return tmp;
}

// Internet version: Fast implmentation of Halton Sequence.
// [[Rcpp::export]]
NumericVector HaltonSeq(const int & k, double & base, int & n) 
{
	NumericVector xk(n);
	int index = k;
	double f;
	for (int i = 0; i < n; i++) {
		index = k + i; 
		f = 1;
		while(index > 0){
			f = f/base;
			xk[i] = xk[i] + f*(fmod(index,base));
			index = index/base;
		}
	}
    return xk;
}

// Fast Implementation of finding x and y order numbers to feed into the linear congruence equation.
// [[Rcpp::export]]
NumericMatrix GetBoxIndices(NumericMatrix& lxy, IntegerVector& base, IntegerVector J)
{
int n = lxy.nrow();
NumericVector ai(2);
NumericMatrix results(n,2);

	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j < 2; j++)
		{
			ai[j] = 0;
			for(int k = 1; k <= J[j]; k++)
			{
				ai[j] += (int(floor(lxy(i,j) * pow( base[j], k))) % base[j]) * pow( base[j], k - 1);
			}
		results(i,j) = ai[j];
		}
	}
	return results;
}

// Given ai values solve the linear congruence.
// [[Rcpp::export]]
NumericVector SolveCongruence(NumericMatrix& A, NumericVector& base, NumericVector J)
{
	NumericVector b = NumericVector::create( pow(base[0], J[0]), pow(base[1], J[1]) );
	double B = b[0]*b[1];

	double j, ll, jj, x;
	bool test;
	NumericVector Index(A.nrow());
	double tmp;
	for(int i = 0; i < A.nrow(); i++)
	{	
		j = 0; test = false;

		ll = (A(i,0) > A(i,1)) ? 0 : 1;		// for speed use the larger value to loop through.
		jj = (ll == 1) ? 0 : 1;				// Track jj as the one to test against.
		while(!test && (j < B))
		{
			x = A(i,ll) + (b[ll])*j;		// Find multiple
			tmp = floor(x / b[jj]);	// Does that particular multiple match with the other?
			test = (x - b[jj]*tmp) == A(i,jj);
			j ++;
		}
		if(test) Index(i) = x;			// Make sure there is an actual solution
		if(!test) Index(i) = -j;		// If no solution make it obvious.
	}
	return Index;
}

// [[Rcpp::export]]
NumericVector SolveCongruence2(NumericMatrix& A, NumericVector& base, NumericMatrix& J)
{
	NumericVector b; 
	double B;

	double j, ll, jj, x;
	bool test;
	NumericVector Index(A.nrow());
	double tmp;
	for(int i = 0; i < A.nrow(); i++)
	{	
		b = NumericVector::create( pow(base[0], J(i,0)), pow(base[1], J(i,1) ));
		B = b[0]*b[1];	
		
		j = 0; test = false;

		ll = (A(i,0) > A(i,1)) ? 0 : 1;		// for speed use the larger value to loop through.
		jj = (ll == 1) ? 0 : 1;				// Track jj as the one to test against.
		while(!test && (j < B))
		{
			x = A(i,ll) + (b[ll])*j;		// Find multiple
			tmp = floor(x / b[jj]);	// Does that particular multiple match with the other?
			test = (x - b[jj]*tmp) == A(i,jj);
			j ++;
		}
		if(test) Index(i) = x;			// Make sure there is an actual solution
		if(!test) Index(i) = -j;		// If no solution make it obvious.
	}
	return Index;
}