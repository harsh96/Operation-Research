
float zero=0.000000001;

double absolute(double a)
{
    if(a<0) return (-1)*a;
    else return a;
}

void swap(double **A,double *B,int i,int j)
{
    float *tmp = A[i];
    A[i]=A[j];
    A[j] = tmp;

    float tmp2 = B[i];
    B[i]=B[j];
    B[j]=tmp2;
}

int isSingular(double **A,double *B,int n,int i)
{
    int j;

    if(absolute(A[i][i])>zero)
    {
        return 0;
    }
    else
    {
        for(j=i+1;j<n;j++)
        {
            if(absolute(A[j][i])>zero)
            {
                swap(A,B,i,j);
                return 0;
            }
        }
        return 1;
    }
}


double * gaussElimination(double **A,double *B,int n)
{
	int i,j,k,flag;
	double factor,*X;

	X = (double *)malloc(n*sizeof(double));


    flag=0;
	for(i=0;i<n-1;i++) //Going From 1st row to (n-1)th row as the base row with base element as i,i
    {
        if(!isSingular(A,B,n,i))
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
        else
        {
            flag=1;
            break;
        }
    	
    }



 	for(i=0;i<n;i++) X[i] = 0;

    if(flag==1)
    {
        return X;
    }

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