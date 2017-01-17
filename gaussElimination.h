double * gaussElimination(double **A,double *B,int n)
{
	int i,j,k;
	double factor,*X;

	X = (double *)malloc(n*sizeof(double));

    //TODO Check if diagonal elements are zero and swap rows

	for(i=0;i<n-1;i++) //Going From 1st row to (n-1)th row as the base row with base element as i,i
    {
    	for(j=i+1;j<n;j++)//Going from the (i+1)th row to the nth row each time with row base element as j,i
    	{
    		factor = A[j][i]/A[i][i];

    		for(k=0;k<n;k++)//For the entire row subtract baserow * (factor)
    		{
    			A[j][k] -= (A[i][k] * factor);
    		}
    		B[j] -= B[i]*factor;
    	}
    }


 	for(i=0;i<n;i++) X[i] = 0;

    for(i=n-1;i>=0;i--)
    {
    	X[i] = B[i];
    	for(j=0;j<i;j++)
    	{
    		X[i] -= A[i][j]*X[j];
    	}
    	X[i] /= A[i][i];
    }

    return X;
}