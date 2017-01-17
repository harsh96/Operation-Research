#include <stdio.h>
#include <stdlib.h>
#include <gaussElimination.h>


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

        flag = 0;
        for(k=0;k<m;k++)
        {
        	if(x[k]!=0) flag=1;
        }
        if(flag!=0)
        {
        	printf("[");
	        for (j=0;j<n;j++)
	        {
	            printf("%.3lf ",solution[j]);
	        }
	        printf("]\n");
        }
        return;
    }

    for (i=start; i<=n-1 && n-i >= m-index; i++)
    {
        combination[index] = i;
        generateSolution(combination,A,B,i+1,index+1,m,n,err);
    }
}


int main()
{
	
	int n,i,j,m,*combination;
	double **A,*B;

	//Get value of N
	printf("Enter the value of N: ");
	scanf("%d",&n);
	printf("Enter the value of M: ");
    scanf("%d",&m);

   	A = (double **)malloc(m *sizeof(double *));
    for (i=0; i<m; i++)
    {
    	A[i] = (double *)malloc(n * sizeof(double));
    }

    B = (double *)malloc(m*sizeof(double));
    combination = (int *)malloc(m*sizeof(int));
    //Taking the MATICES as input
    printf("Enter MATRIX A\n");
    for(i=0;i<m;i++)
    {
    	for(j=0;j<n;j++)
    	{
    		scanf("%lf",&A[i][j]);
    	}
    }

    printf("Enter MATRIX B\n");
    for(i=0;i<n;i++)
    {
		scanf("%lf",&B[i]);
    }

    printf("Printing Basic solutions \n");
    generateSolution(combination,A,B,0,0,m,n);

    X = gaussElimination(A,B,n);


    for(i=0;i<n;i++) printf("%lf ",X[i]);
    	

	return 0;
}