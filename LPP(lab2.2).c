#include <stdio.h>
#include <stdlib.h>

float zero=0.000000001; //as floating point equality 0 doesn’t work
double **data;
int count=0;


double absolute(double a)
{
    if(a<0) return (-1)*a;
    else return a;
}

int NCM(int n,int m)
{
    int temp = 1,i;
    for(i=n;i>(n-m);i--)
    {
        temp *= i;
    }
    for(i=m;i>1;i--)
    {
        temp /= i;
    }
    return temp;
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


void generateSolution(int combination[],double **A,double *B,int start,int index,int m,int n)
{
    int j=0,k=0,i,flag;
    double *x,*solution;
    x=(double *)malloc(m*sizeof(double));
    solution=(double *)malloc(n*sizeof(double));

    if(index==m)
    {
        double **A_square;
        A_square=(double ** )malloc(m*sizeof(double * ));
        for (j=0;j<m;j++)
        {
            A_square[j]=(double *)malloc(m*sizeof(double));
        }
        for (j=0;j<m;j++)
        {
            for (k=0;k<m;k++)
            {

                A_square[j][k]=A[j][combination[k]];
            }
        }

        x = gaussElimination(A_square,B,m);

        for(j=0;j<n;j++)
        {
            solution[j]=0;
        }
        for (k=0;k<m;k++)
        {
            solution[combination[k]]=x[k];
        }

        data[count] = solution;
        count++;

        return;
    }

    for (i=start; i<=n-1 && n-i >= m-index; i++)
    {
        combination[index] = i;
        generateSolution(combination,A,B,i+1,index+1,m,n);
    }
}


int main()
{

	int n,i,j,m,*combination,input,flag,temp;
	double **A,*B,*expression;

	//Get value of N
	printf("Enter the number of variables: ");
	scanf("%d",&n);
	printf("Enter the number of equations: ");
    scanf("%d",&m);

    data = (double **)malloc(NCM(n+m,m) * sizeof(double *));
    for (i=0; i<n; i++)
    {
        data[i] = (double *)malloc((n+m) * sizeof(double));
    }

   	A = (double **)malloc(m * sizeof(double *));
    for (i=0; i<m; i++)
    {
    	A[i] = (double *)malloc((n+m) * sizeof(double));
    }

    B = (double *)malloc(m*sizeof(double));
    expression = (double *)malloc(n*sizeof(double));
    combination = (int *)malloc(m*sizeof(int));
    //Taking the MATICES as input
    printf("Enter LHS Matrix \n");
    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            A[i][j] =1;
        }
    }
    for(i=0;i<m;i++)
    {
    	for(j=0;j<n;j++)
    	{
    		scanf("%lf",&A[i][j]);
    	}
    }

    printf("Enter RHS Matrix\n");
    for(i=0;i<m;i++)
    {
		scanf("%lf",&B[i]);
    }

    generateSolution(combination,A,B,0,0,m,n+m);

    printf("Enter the coefficients of expression to be maximised\n");
    for(i=0;i<n;i++)
    {
        scanf("%lf",&expression[i]);
    }
    

	return 0;
}
