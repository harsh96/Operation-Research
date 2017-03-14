#include <stdio.h>
#include <stdlib.h>
#include <string.h>

float zero=0.000000001; //as floating point equality 0 doesn’t work
int unbounded = 0;
int infeasible = 0;
int alternate = 0;
int altIdx = 0;

double absolute(double a)
{
    if(a<0) return (-1)*a;
    else return a;
}

void makeTableau(double **A,double *B,double *expression,int rows,int columns,double **tableau,int **additional,int n)
{
    int i,j,k,m,dummy1,dummy2;
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

        for(j=0;j<m;j++)
        {
            if(additional[1][j]==1)
            {
                tableau[1][i+2] += (-1)*A[j][i];
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
        tableau[1][n+i+2] = 1;  //z-row
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

    dummy1 = 0;
    dummy2 = 0;

    for(i=(surplus+n+2);i<columns-1;i++)
    {
        k= i-(surplus+n+2);

        if(additional[0][k] == 1)
        {
            dummy1++;
            tableau[0][i] = 200+dummy1;//slack
            tableau[k+2][i] = 1;
            tableau[k+2][0] = 200+dummy1;
        }
        else if(additional[1][k] == 1)
        {
            dummy2++;
            tableau[0][i] = 300+dummy2;//artificial
            tableau[k+2][i] = 1 ;
            tableau[k+2][0] = 300+dummy2;
        }
    }

    for(i=0;i<m;i++)
    {
        if(additional[1][i]==1) tableau[1][columns-1] -= B[i];
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
            printf("Iteration %d:\n",totalIterations);
            performSimplexIteration(tableau,rows,columns,colIdx);
            printTableau(tableau,rows,columns);
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

void makeTableauPhase2(double **tableau,double **tableauPhase2,double *expression,int rows,int columns,int n)
{
    //here rows and columns are number of rows and columns of tableau
    int i,j,k,m,temp;
    double factor;

    for(i=0;i<rows;i++)
    {
        tableauPhase2[i][0] = tableau[i][0];
        tableauPhase2[i][1] = tableau[i][1];
    }
    k=2;
    for(i=2;i<columns;i++)
    {
        if(tableau[0][i]<300)
        {
            for(j=0;j<rows;j++)
            {
                tableauPhase2[j][k] = tableau[j][i];
                tableauPhase2[1][k] = 0;
            }
            k++;
        }
    }

    for(i=0;i<n;i++)
    {
        tableauPhase2[1][i+2] = (-1)*expression[i]*expression[n];
    }

    for(i=1;i<=n;i++)
    {
        factor = tableauPhase2[1][i+1];
        for(j=2;j<rows;j++)
        {
            temp = tableauPhase2[j][0];
            if(temp==i)
            {
                for(m=2;m<k;m++)
                {
                    tableauPhase2[1][m] -= factor*tableauPhase2[j][m];
                }
            }
        }
    }
}

int main()
{

	int n,i,j,m,input,flag,temp,totalVariables,artificial;
    char symbol,dummy;
	double **A,*B,*expression,**tableau,**tableauPhase2;
    int **additional;
	printf("Simplex Method\n");
	printf("------------------------------\n");
	//Get value of N
	printf("Enter the number of variables: ");
	scanf("%d",&n);
	printf("Enter the number of constraint equations: ");
    scanf("%d",&m);

    printf("Enter Constraint equations with < or > or =\n");

    A = (double **)malloc(m * sizeof(double *));
    for (i=0; i<m; i++)
    {
        A[i] = (double *)malloc(n * sizeof(double));
    }

    B = (double *)malloc(m*sizeof(double));

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

        scanf("%lf",&B[i]);

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

    expression = (double *)malloc((n+1)*sizeof(double));


    printf("Enter the coefficients of expression and 1 if you want to maximise else -1\n");
    for(i=0;i<=n;i++)
    {
        scanf("%lf",&expression[i]);
    }

    printf("PHASE I:\n");
    makeTableau(A,B,expression,m+2,totalVariables+3,tableau,additional,n);
    printf("Initial tableau:\n");
    printTableau(tableau,m+2,totalVariables+3);
    applySimplex(tableau,m+2,totalVariables+3);
    printf("Final tableau:\n");
    printTableau(tableau,m+2,totalVariables+3);

    if(absolute(tableau[1][totalVariables+2])>zero) infeasible = 1;
    else
    {
        artificial = 0;
        for(i=0;i<m;i++)
        {
            if(additional[1][i]==1) artificial++;
        }
        tableauPhase2 = (double **)malloc((m+2)*sizeof(double *));
        for(i=0; i<m+2; i++)
        {
            tableauPhase2[i] = (double *)malloc((totalVariables+3-artificial) * sizeof(double));
        }
        printf("PHASE II:\n");
        makeTableauPhase2(tableau,tableauPhase2,expression,m+2,totalVariables+3,n);
        printf("Initial tableau:\n");
        printTableau(tableauPhase2,m+2,totalVariables+3-artificial);
        applySimplex(tableauPhase2,m+2,totalVariables+3-artificial);
        printf("Final tableau:\n");
        printTableau(tableauPhase2,m+2,totalVariables+3-artificial);
    }
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
        printf("%.5lf\n",(tableauPhase2[1][totalVariables+2-artificial])*expression[n]);
        printf("for values of xi as\n");
        printf("[ ");
        for(i=1;i<=n;i++)
        {
            for(j=2;j<m+2;j++)
            {
                temp = tableauPhase2[j][0];//typecasting to int
                if(i == temp)
                {
                    printf("%.3lf ",tableauPhase2[j][totalVariables+2-artificial]);
                    break;
                }
            }
            if(j==(m+2)) printf("0.000 ");

        }
        printf("]\n");
        if(alternate==1)
        {
            printf("Alternate solution also exists and it will be the convex of above and \n");
            performSimplexIteration(tableauPhase2,m+2,totalVariables+3-artificial,altIdx);
            printf("[ ");
            for(i=1;i<=n;i++)
            {
                for(j=2;j<m+2;j++)
                {
                    temp = tableauPhase2[j][0];//typecasting to int
                    if(i == temp)
                    {
                        printf("%.3lf ",tableauPhase2[j][totalVariables+2-artificial]);
                        break;
                    }
                }
                if(j==(m+2)) printf("0.000 ");
            }
            printf("]\n");
        }
    }

	return 0;
}
