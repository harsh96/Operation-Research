#include <stdio.h>
#include <stdlib.h>
#include <string.h>

float zero=0.000000001; //as floating point equality 0 doesn’t work
int unbounded = 0;
int infeasible = 0;
int alternate = 0;
int altIdx = 0;
double M = 10000;

double absolute(double a)
{
    if(a<0) return (-1)*a;
    else return a;
}

void makeTableau(double **A,double *B,double *expression,int n,int m,double **tableau)
{
    int i,j,k;
    int rows = m+2, columns = m+n+2;
    for(i=0;i<rows;i++)
    {
        for(j=0;j<columns;j++)
        {
            tableau[i][j] = 0;
        }
    }

    for(i=1;i<=n;i++)
    {
        tableau[0][i] = i; //variable number
        // For the z-row
        tableau[1][i] = (-1)*expression[i-1]*expression[n];//expression[n] is 1 if minimisation
        //z-row ends

        for(j=0;j<m;j++)
        {
            tableau[j+2][i] = A[j][i-1];
        }
    }

    for(i=n+1;i<=n+m;i++)
    {
        tableau[i-n+1][i] = 1;
        tableau[0][i] = 200+i-n;
        tableau[i-n+1][0] = 200+i-n;
    }

    for(i=2;i<rows;i++)
    {
        tableau[i][columns-1] = B[i-2];
    }
}

void printTableau(double **tableau,int rows, int columns)
{
    int i,j,temp;
    printf("  Basic  ");
    for(j=1;j<columns-1;j++)
    {
        temp = ((int)tableau[0][j])%100;
        if(tableau[0][j] > 200)
        {
            printf("  slack%d ",temp);
        }
        else
        {
            printf("    X%d   ",temp);
        }
    }
    printf("Solution \n");
    for(i=1;i<rows;i++)
    {
        if(i==1) 
        {
            printf("    Z    ");
        }
        else
        {
            temp = ((int)tableau[i][0])%100;
            if(tableau[i][0] > 200)
            {
                printf("  slack%d ",temp);
            }
            else
            {
                printf("    X%d   ",temp);
            }
        }
        for(j=1;j<columns;j++)
        {   
            printf("%8.2lf ",tableau[i][j]);
        }
        printf("\n");
    }
}

int performDualSimplexIteration(double **tableau,int rows,int columns)
{
    int i,j,rowIdx,colIdx;
    double minRatio,pivot,entCol,minNeg;

    minNeg = 0;
    minRatio = 1000000;

    for(i=2;i<rows;i++)
    {
        if(tableau[i][columns-1]<minNeg)
        {
            minNeg = tableau[i][columns-1];
            rowIdx = i;
        }
    }
    if(minNeg==0) return 0; //Optimality reached
    else
    {
        for(j=1;j<columns-1;j++)
        {
            if(tableau[rowIdx][j]<0)
            {
                if((tableau[1][j]/tableau[rowIdx][j])<minRatio)
                {
                    minRatio = tableau[1][j]/tableau[rowIdx][j];
                    colIdx = j;
                }
            }
        }
        if(minRatio==1000000) return -1; //Not feasible
        else
        {
            pivot = tableau[rowIdx][colIdx]; 

            //new pivot equation
            for(j=1;j<columns;j++)
            {
                tableau[rowIdx][j] = tableau[rowIdx][j]/pivot;
            }

            //other new equations

            for(i=1;i<rows;i++)
            {
                if(i!=rowIdx)
                {
                    entCol = tableau[i][colIdx];
                    for(j=1;j<columns;j++)
                    {
                        tableau[i][j] -= entCol * tableau[rowIdx][j];
                    }
                }
            }

            tableau[rowIdx][0] = tableau[0][colIdx];
            return 1;// iteration complete
        }
    }

}

int applyDualSimplex(double **tableau,int rows,int columns)
{
    int temp,totalIterations;
    totalIterations = 0;
    while(1)
    {
        temp = performDualSimplexIteration(tableau,rows,columns);
        if(temp!=1) break;
        else totalIterations++;
    }
    if(temp==-1) infeasible = 1;

    return totalIterations;
}

int main()
{

	int n,i,j,m,input,temp,totalIterations;
    char symbol,dummy;
	double **A,*B,*expression,**tableau,**dummyTableau;
	printf("Dual Simplex Method\n");
	printf("------------------------------\n");
	printf("Enter the number of variables: ");
	scanf("%d",&n);
	printf("Enter the number of constraint equations: ");
    scanf("%d",&m);

    printf("Enter Constrain Equation followed by < or > \n");

    A = (double **)malloc(m * sizeof(double *));
    for(i=0; i<m; i++)
    {
        A[i] = (double *)malloc(n * sizeof(double));
    }

    B = (double *)malloc(m*sizeof(double));

    for(i=0;i<m;i++)
    {
    	for(j=0;j<n;j++)
    	{
    		scanf("%lf",&A[i][j]);
    	}

        scanf("%c",&dummy);
        scanf("%c",&symbol);

        scanf("%lf",&B[i]);

        if(symbol == '>')
        {
            for(j=0;j<n;j++)
            {
                A[i][j] *= -1;
            } 
            B[i] *= -1;
        }
    }

    tableau = (double **)malloc((m+2) * sizeof(double *));
    for (i=0; i<m+2; i++)
    {
        tableau[i] = (double *)malloc((n+m+2) * sizeof(double));
    }

    dummyTableau = (double **)malloc((m+2) * sizeof(double *));
    for (i=0; i<m+2; i++)
    {
        dummyTableau[i] = (double *)malloc((n+m+2) * sizeof(double));
    }
    
    expression = (double *)malloc((n+1)*sizeof(double));

    printf("Enter the coefficients of expression and 1 if you want to maximise else -1\n");
    for(i=0;i<=n;i++)
    {
        scanf("%lf",&expression[i]);
    }
    expression[n] *= -1 ;//since we need 1 for minimisation

    makeTableau(A,B,expression,n,m,tableau);
    totalIterations = applyDualSimplex(tableau,m+2,m+n+2);

    printf("----------------------------\n");
    printf("Enter 1 to view intitial tableau\n");
    printf("Enter 2 to print table of ith iteration\n");
    printf("Enter 3 to view the optimal solution\n");
    printf("Enter 4 to view the final tableau\n");
    printf("Enter 5 to Exit\n");
    while(1)
    {
        printf("----------------------------\n");
        scanf("%d",&input);
        if(input == 1)
        {
            makeTableau(A,B,expression,n,m,dummyTableau);
            printTableau(dummyTableau,m+2,m+n+2);
        }
        else if(input==2)
        {
            printf("Enter the value of i: ");
            scanf("%d",&i);
            if(i<totalIterations)
            {
                makeTableau(A,B,expression,n,m,dummyTableau);
                while(i--)
                {
                   performDualSimplexIteration(dummyTableau,m+2,m+n+2); 
                }
                printTableau(dummyTableau,m+2,m+n+2);
            }
            else
            {
                printf("i cannot be greater than %d\n",totalIterations);
            }
        }
        else if(input==3)
        {
            if(unbounded==1)
            {
                printf("The expression is unbounded\n");
            }
            else if(infeasible==1)
            {
                printf("The LPP is infeasible\n");
            }
            else 
            {
                
                if(expression[n]==-1)
                {
                    printf("The maximum value of the expression "); 
                }
                else
                {
                    printf("The minimum value of the expression ");
                }
                for(i=0;i<n-1;i++)
                {
                    printf("(%.2lf)x%d + ",expression[i],(i+1));
                }
                printf("(%.2lf)x%d is ",expression[n-1],(n));
                printf("%.5lf\n",(tableau[1][m+n+1])*expression[n]);
                printf("for values of xi as\n");
                printf("[ ");
                for(i=1;i<=n;i++)
                {
                    for(j=2;j<m+2;j++)
                    {
                        temp = tableau[j][0];//typecasting to int
                        if(i == temp)
                        {
                            printf("%.3lf ",tableau[j][m+n+1]);
                            break;
                        }
                    }
                    if(j==(m+2)) printf("0.000 ");

                }
                printf("]\n");
                
                
            }
        }
        else if(input==4)
        {
            printTableau(tableau,m+2,m+n+2);
        }
        else break;
    }

    
	return 0;
}
