#include <stdio.h>
#include <stdlib.h>

float zero=0.000000001; //as floating point equality 0 doesn’t work
int unbounded = 0;

double absolute(double a)
{
    if(a<0) return (-1)*a;
    else return a;
}

void makeTableau(double **A,double *B,double *expression,int rows,int columns,double **tableau)
{
    int n,i,j;
    n = columns-rows-1;//1 to n are given variables more than these are additional

    for(i=0;i<rows;i++)
    {
        for(j=0;j<columns;j++)
        {
            tableau[i][j] = 0;
        }
    }

    for(j=0;j<n;j++)
    {
        tableau[1][2+j] = expression[j];
    }

    for(i=2;i<rows;i++)
    {
        tableau[i][0] = n+(i-1); //starts from n+1 denoting the variable number
    }
    for(i=2;i<columns-1;i++)
    {
        tableau[0][i] = i-1; //all variable numbers 1 to n+m(since only slack)
    }
    for(i=2;i<rows;i++)
    {
        for(j=2;j<columns-1;j++)
        {

            tableau[i][j] = A[i-2][j-2];
        }
        tableau[i][columns-1] = B[i-2];
    }
    printf("\n");
}

void printTableau(double **tableau,int rows, int columns)
{
    int i,j;
    for(i=0;i<rows;i++)
    {
        for(j=0;j<columns;j++)
        {
            printf("%.2lf ",tableau[i][j]);
        }
        printf("\n");
    }
}

void applySimplex(double **tableau,int rows,int columns)
{
    int i,j,colIdx,rowIdx;
    double minNeg,minRatio,pivot,entCol;
    while(1)
    {
        minNeg = 0;
        for(i=2;i<columns-1;i++)
        {
            if(tableau[1][i]<minNeg)
            {
                minNeg = tableau[1][i];
                colIdx = i;
            }
        }

        if(absolute(minNeg)<zero) break;//all positive reached optimal state
        for(i=2;i<rows;i++)
        {
            if(tableau[i][colIdx] > 0) break;
        }
        if(i==rows)
        {
            //unbounded
            unbounded = 1;
            break;
        }
        else
        {
            minRatio = 1000000;
            for(i=2;i<rows;i++)
            {
                if(tableau[i][colIdx] > 0) 
                {
                    if((tableau[i][columns-1]/tableau[i][colIdx])<minRatio)
                    {
                        minRatio = tableau[i][columns-1]/tableau[i][colIdx];
                        rowIdx = i;
                    }
                }
            }

            pivot = tableau[rowIdx][colIdx];

            //new pivot equation
            for(i=2;i<columns;i++)
            {
                tableau[rowIdx][i] = tableau[rowIdx][i]/pivot;
            }

            //other new equations

            for(i=1;i<rows;i++)
            {
                if(i!=rowIdx)
                {
                    entCol = tableau[i][colIdx];
                    for(j=2;j<columns;j++)
                    {
                        tableau[i][j] -= entCol * tableau[rowIdx][j];
                    }
                }
            }

            tableau[rowIdx][0] = tableau[0][colIdx];
        }
    }
}

int main()
{

	int n,i,j,m,input,flag,temp;
	double **A,*B,*expression,**tableau;
	printf("This program will maximize an expression with constraints of the form <=\n");
	printf("------------------------------\n");
	//Get value of N
	printf("Enter the number of variables: ");
	scanf("%d",&n);
	printf("Enter the number of constraint equations: ");
    scanf("%d",&m);


   	A = (double **)malloc(m * sizeof(double *));
    for (i=0; i<m; i++)
    {
    	A[i] = (double *)malloc((n+m) * sizeof(double));
    }

    tableau = (double **)malloc((m+2) * sizeof(double *));
    for (i=0; i<m+2; i++)
    {
        tableau[i] = (double *)malloc((n+m+3) * sizeof(double));
    }

    B = (double *)malloc(m*sizeof(double));
    expression = (double *)malloc(n*sizeof(double));

    printf("Enter LHS Constraint Matrix \n");
    
    for(i=0;i<m;i++)
    {
    	for(j=0;j<n;j++)
    	{
    		scanf("%lf",&A[i][j]);
    	}
    }

    //adding slack variables
    for(i=0;i<m;i++)
    {
        for(j=n;j<n+m;j++)
        {
            A[i][j] =0;
        }
        A[i][i+n] =1;
    }

    printf("Enter RHS Constraint Matrix\n");
    for(i=0;i<m;i++)
    {
		scanf("%lf",&B[i]);
    }

    printf("Enter the coefficients of expression to be maximised\n");
    for(i=0;i<n;i++)
    {
        scanf("%lf",&expression[i]);
        expression[i] = -1*expression[i];
    }

    makeTableau(A,B,expression,m+2,n+m+3,tableau);
    applySimplex(tableau,m+2,n+m+3);
    if(unbounded == 0)
    {
        printf("The maximum value of the expression ");
        for(i=0;i<n-1;i++)
        {
            printf("(%.2lf)x%d + ",(-1)*expression[i],(i+1));
        }
        printf("(%.2lf)x%d is ",(-1)*expression[n-1],(n));
        
        printf("%.5lf\n",tableau[1][m+n+2]);
        printf("for values of xi as\n");
        printf("[ ");
        for(i=1;i<=n;i++)
        {
            for(j=2;j<m+2;j++)
            {
                temp = tableau[j][0];//typecasting to int
                if(i == temp)
                {
                    printf("%.3lf ",tableau[j][m+n+2]);
                    break;
                }
            }
        }
        printf("]\n");  
    }
    else
    {
        printf("This set is unbounded\n");
    }
    
    
	return 0;
}
