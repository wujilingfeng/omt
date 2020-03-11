#ifndef Cholesky_ALG
#define Cholesky_ALG
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
void cholesky_decomp(double** a,double **L,double *D,int n)
{
    for(int k=0;k<n;k++)
	{
		double sum = 0;
        //j rows
		for (int j = k;j<n;j++)
		{
			sum = 0;
			for (int m=0;m<k;m++)
			{
				if (m==k)
				{
					continue;
				}
				else
				{
					sum+= L[j][m] * D[m] * L[k][m];
				}
			}
			if (fabs(a[j][k] - sum) <= 1e-10)
			{
				if (j == k) {
					D[k] = 0;
				}
				else {
					L[j][k] = 0;
				}
			}
			else
			{
				if (j == k) {

					D[k] = a[j][k] - sum;
				}
				else {

					L[j][k] = (a[j][k] - sum) / double(D[k]);
				}
			}
		}
	}
}
double *cholesky_solve(double** a, double** L, double* D, double* b, int n)
{
	double *y =(double*)malloc(sizeof(double)*n);
	double *x =(double*)malloc(sizeof(double)*n);
	double sum = 0;
	for (int i = 0;i<n;i++)
	{
		sum = 0;
		int k = 0;
		while (k<i)
		{
			sum += D[k] * L[i][k] * y[k];
			k++;
		}
		y[i] = (b[i] - sum) / double(L[i][i] * D[i]);

	}
	for (int i = n - 1;i >= 0;i--)
	{
		sum = 0;
		int k = n - 1;
		while (k>i)
		{
			sum += L[k][i] * x[k];
			k--;
		}
		//assert(fabs(L[i][i]) >= 1e-6);
		x[i] = (y[i] - sum) / double(L[i][i]);
	}
    free(y);
	return x;
}
#endif

