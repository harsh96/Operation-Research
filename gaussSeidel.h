
int isDiagonallyDominant(double **A,int n)
{
    int i,j;
    int flag=1;
    double sum=0,temp,temp2;
    for (i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {

            if(j!=i){

                if(A[i][j]<0) temp = (-1)*A[i][j];
                else temp = A[i][j];

                sum += temp;
            }
        }

        if(A[i][i]<0) temp2 = (-1)*A[i][i];
        else temp2 = A[i][i];

        if(temp2 < sum)
        {
            flag=0;
            return flag;
        }
    }
    return flag;
}


double * gaussSeidel(double **A,double *B,int n,float err)
{
    int key=1,i,j;
    double *X,*X1,sum;
    float calErr;

    X = (double *)malloc(n*sizeof(double));
    X1 = (double *)malloc(n*sizeof(double));

    for(i=0;i<n;i++)
    {
        X[i] = 0;
    }

    //TODO attempt to make it diagonally dominant.

    if (!isDiagonallyDominant(A,n))
    {
        key=0;
        printf("Not diagonally dominant, method fails->");
    }

    while(key==1)
    {

        for(i=0;i<n;i++)
        {
            sum = B[i];
            for(j=0;j<n;j++)
            {
                if(i!=j)
                {
                  sum -= A[i][j]*X[j];
                }
            }
            X1[i] = sum/A[i][i];
            calErr = (X1[i]-X[i])/X1[i];
            if(calErr<0) calErr *= (-1);
            if(calErr>err)
            {
                key = 1;
                X[i] = X1[i];
            }
            else
            {
                key = 0;
            }
        }
    }
    return X1;
}