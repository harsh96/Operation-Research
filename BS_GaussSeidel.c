#include <stdio.h>
#include <stdlib.h>
#include <gaussSeidel.h>

void generateSolution(int combination[],double **A,double *B,int start,int index,int m,int n,double err)
{
    int j=0,k=0,i;
    double *x,*solution;
    x=(double *)malloc(m*sizeof(double));
    solution=(double *)malloc(n*sizeof(double));

    if (index==m)
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

        x = gaussSeidel(A_square,B,m,err);

        for(j=0;j<n;j++)
        {
            solution[j]=0;
        }
        for (k=0;k<m;k++)
        {
            solution[combination[k]]=x[k];
        }
        printf("[");
        for (j=0;j<n;j++)
        {
            printf("%.3lf ",solution[j]);
        }
        printf("]\n");
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

	int n,m,i,j,key,*combination;
    float err=0.001;
	double **A,*B;

	printf("Enter the value of N: ");
	scanf("%d",&n);
    printf("Enter the value of M: ");
    scanf("%d",&m);
    printf("Enter the value of error permissible: ");
    scanf("%f",&err);

   	A = (double **)malloc(m *sizeof(double *));
    for (i=0; i<m; i++)
    {
    	A[i] = (double *)malloc(n * sizeof(double));
    }

    B = (double *)malloc(m*sizeof(double));

    combination = (int *)malloc(m*sizeof(int));

    printf("Enter MATRIX A\n");
    for(i=0;i<m;i++)
    {
    	for(j=0;j<n;j++)
    	{
    		printf("A[%d][%d]=",i,j);
    		scanf("%lf",&A[i][j]);
    	}
    }

    printf("Enter MATRIX B\n");
    for(i=0;i<m;i++)
    {
    	printf("B[%d]=",i);
		scanf("%lf",&B[i]);
    }
    printf("Printing Basic solutions \n");
    generateSolution(combination,A,B,0,0,m,n,err);

	return 0;
}
