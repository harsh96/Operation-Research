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

void makeTableau(double **A,double *B,double *expression,int rows,int columns,double **tableau,int **additional,int n)
{
    int i,j,k,m;
    int slack=0,surplus=0,artificial=0;

    m = rows-2;

    for(i=0;i<m;i++)
    {
        if(additional[0][i]==1) slack++;
        if(additional[0][i]==-1) surplus++;
        if(additional[1][i]==1) artificial++;
    }

    for(i=0;i<rows;i++)
    {
        for(j=0;j<columns;j++)
        {
            tableau[i][j] = 0;
        }
    }
    //Order in tableau is
    // Normal | Surplus | Slack and Artificial
    // 1-n | 100+ | 200+ and 300+

    for(i=0;i<n;i++)
    {
        tableau[0][i+2] = i+1; //variable number
        // For the z-row
        tableau[1][i+2] = (-1)*expression[i]*expression[n];
        for(j=0;j<m;j++)
        {
            if(additional[1][j]==1)
            {
                tableau[1][i+2] += M*(-1)*A[j][i]; 
            }
        }
        //z-row ends

        for(j=0;j<m;j++)
        {
            tableau[j+2][i+2] = A[j][i];
        }
    }

    for(i=0;i<surplus;i++)
    {
        tableau[0][n+i+2] = 100+i+1; //variable number
        tableau[1][n+i+2] = M;  //z-row
        k=0;
        for(j=0;j<m;j++)
        {
            if(additional[0][j]==-1) k++;
            if(k==(i+1))
            {
                tableau[j+2][n+i+2] = -1;
                break;
            }
            //else it will be zero which was initialised
        }
    }

    for(i=(surplus+n+2);i<columns-1;i++)
    {
        k= i-(surplus+n+2);

        if(additional[0][k] == 1)
        {
            tableau[0][i] = 200+k+1;//slack
            tableau[k+2][i] = 1;
            tableau[k+2][0] = 200+k+1;
        }
        else if(additional[1][k] == 1)
        {
            tableau[0][i] = 300+k+1;//artificial
            tableau[k+2][i] = 1 ;
            tableau[k+2][0] = 300+k+1;
        }
    }

    for(i=0;i<m;i++)
    {
        if(additional[1][i]==1) tableau[1][columns-1] -= M*B[i];
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
    for(j=2;j<columns-1;j++)
    {
        temp = ((int)tableau[0][j])%100;
        if(tableau[0][j] > 300)
        {
           printf("  artif%d ",temp); 
        }
        else if(tableau[0][j] > 200)
        {
            printf("  slack%d ",temp);
        }
        else if(tableau[0][j]>100)
        {
            printf(" surplus%d ",temp);
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
            if(tableau[i][0] > 300)
            {
               printf("  artif%d ",temp); 
            }
            else if(tableau[i][0] > 200)
            {
                printf("  slack%d ",temp);
            }
            else if(tableau[i][0]>100)
            {
                printf(" surplus%d ",temp);
            }
            else
            {
                printf("    X%d   ",temp);
            }
        }
        for(j=2;j<columns;j++)
        {   
            printf("%8.2lf ",tableau[i][j]);
        }
        printf("\n");
    }
}

void performSimplexIteration(double **tableau,int rows,int columns,int colIdx)
{
    int i,j,rowIdx;
    double minRatio,pivot,entCol;

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

int applySimplex(double **tableau,int rows,int columns)
{
    int i,j,colIdx,totalIterations=0;
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
            totalIterations++;
            performSimplexIteration(tableau,rows,columns,colIdx);
        }    
    }
    for(i=2;i<rows;i++)
    {
        if(tableau[i][0]>300)//artificial
        {
            if(absolute(tableau[i][columns-1])>zero)
            {
                infeasible = 1;
            }
        }
    }
    if(infeasible != 1)
    {
        for(i=2;i<columns-1;i++)
        {
            if(absolute(tableau[1][i])<zero)
            {
                for(j=2;j<rows;j++)
                {
                    if(tableau[0][i]==tableau[j][0]) break;
                }
                if(j==rows)
                {
                    alternate = 1;
                    altIdx = i;
                }
            }
        }
    }
    return totalIterations;
}

int main()
{

	int n,i,j,m,input,flag,temp,totalVariables;
    char symbol,dummy;
	double **A,*B,*expression,**tableau;
    int **additional;
	printf("Simplex Method\n");
	printf("------------------------------\n");
	//Get value of N
	printf("Enter the number of variables: ");
	scanf("%d",&n);
	printf("Enter the number of constraint equations: ");
    scanf("%d",&m);

    printf("Enter LHS Constraint Matrix followed by < or > or =\n");

    A = (double **)malloc(m * sizeof(double *));
    for (i=0; i<m; i++)
    {
        A[i] = (double *)malloc(n * sizeof(double));
    }

    additional = (int **)malloc(2 * sizeof(int *));
    additional[0] = (int *)malloc(m * sizeof(int)); //Stores slack or surplus variable
    additional[1] = (int *)malloc(m * sizeof(int)); //Stores artificial variable
    totalVariables = n;

    for(i=0;i<m;i++)
    {
    	for(j=0;j<n;j++)
    	{
    		scanf("%lf",&A[i][j]);
    	}

        scanf("%c",&dummy);
        scanf("%c",&symbol);

        if(symbol == '<')
        {
            additional[0][i] = 1;
            additional[1][i] = 0;
            totalVariables++;
        }
        else if(symbol == '>')
        {
            additional[0][i] = -1;
            additional[1][i] = 1;
            totalVariables++;
            totalVariables++;
        }
        else
        {
            additional[0][i] = 0;
            additional[1][i] = 1;
            totalVariables++;
        }
    }

    tableau = (double **)malloc((m+2) * sizeof(double *));
    for (i=0; i<m+2; i++)
    {
        tableau[i] = (double *)malloc((totalVariables+3) * sizeof(double));
    }

    B = (double *)malloc(m*sizeof(double));
    expression = (double *)malloc((n+1)*sizeof(double));

    printf("Enter RHS Constraint Matrix\n");
    for(i=0;i<m;i++)
    {
		scanf("%lf",&B[i]);
    }

    printf("Enter the coefficients of expression and 1 if you want to maximise else -1\n");
    for(i=0;i<=n;i++)
    {
        scanf("%lf",&expression[i]);
    }

    makeTableau(A,B,expression,m+2,totalVariables+3,tableau,additional,n);
    applySimplex(tableau,m+2,totalVariables+3);
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
        
        if(expression[n]==1)
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
        printf("%.5lf\n",(tableau[1][totalVariables+2])*expression[n]);
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
            if(j==(m+2)) printf("0.000 ");

        }
        printf("]\n");
        if(alternate==1)
        {
            printf("Alternate solution also exists and it will be the convex of above and \n");
            performSimplexIteration(tableau,m+2,totalVariables+3,altIdx);
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
        printTableau(tableau,m+2,totalVariables+3);
    }

	return 0;
}
