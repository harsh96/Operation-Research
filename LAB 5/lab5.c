#include <stdio.h>
#include <stdlib.h>
#include <string.h>

float zero=0.000000001; //as floating point equality 0 doesn’t work
int unbounded = 0;
int infeasible = 0;
int alternate = 0;
int altIdx = 0;
double M = 10000;
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
    double *tmp = A[i];
    A[i]=A[j];
    A[j] = tmp;

    double tmp2 = B[i];
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
        for(j=n-1;j>i;j--)
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
        double **A_square,*b_temp;
        b_temp = (double *)malloc(m*sizeof(double));
        A_square=(double ** )malloc(m*sizeof(double * ));
        for (j=0;j<m;j++)
        {
            A_square[j]=(double *)malloc(m*sizeof(double));
            b_temp[j] = B[j];
        }
        for (j=0;j<m;j++)
        {
            for (k=0;k<m;k++)
            {

                A_square[j][k]=A[j][combination[k]];
            }
        }


        x = gaussElimination(A_square,b_temp,m);

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

void makeTableau(double **A,double *B,double *expression,int rows,int columns,double **tableau,int **additional,int n)
{
    int i,j,k,m;
    int dummy1, dummy2;
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
        if(tableau[0][j] > 300) printf("  artif%d ",temp);
        else if(tableau[0][j] > 200) printf("  slack%d ",temp);
        else if(tableau[0][j]>100) printf(" surplus%d ",temp);
        else printf("    X%d   ",temp);
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
            if(tableau[i][0] > 300) printf("  artif%d ",temp); 
            else if(tableau[i][0] > 200) printf("  slack%d ",temp);
            else if(tableau[i][0]>100) printf(" surplus%d ",temp); 
            else printf("    X%d   ",temp);
        }
        for(j=2;j<columns;j++)
        {   
            printf("%8.2lf ",tableau[i][j]);
        }
        printf("\n");
    }
}

int performSimplexIteration(double **tableau,int rows,int columns,int idx)
{
    int i,j,colIdx,rowIdx;
    double minNeg,minRatio,pivot,entCol;
    
    minNeg = 0;
    if(idx != -1)
    {
        colIdx = idx;
    }
    else
    {
        for(i=2;i<columns-1;i++)
        {
            if(tableau[1][i]<minNeg)
            {
                minNeg = tableau[1][i];
                colIdx = i;
            }
        }
        if(absolute(minNeg)<zero) return 0;//all positive reached optimal state
    }
    
    for(i=2;i<rows;i++)
    {
        if(tableau[i][colIdx] > 0) break;
    }
    if(i==rows)
    {
        //unbounded
        return -1;
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
        return 1;// iteration complete
    }    


    
}

int applySimplex(double **tableau,int rows,int columns)
{
    int totalIterations=0,temp,i,j;

    while(1)
    {
        temp = performSimplexIteration(tableau,rows,columns,-1);
        if(temp!=1) break;
        else totalIterations++;       
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

	int n,i,j,m,input,flag,temp,totalVariables,*combination,totalIterations;
    char symbol,dummy;
	double **A,*B,*expression,**tableau,**dummyTableau,**A_BFS;
    int **additional;
	printf("Operation Research Lab 5\n");
	printf("------------------------------\n");
	//Get value of N
	printf("Enter the number of variables: ");
	scanf("%d",&n);
	printf("Enter the number of constraint equations: ");
    scanf("%d",&m);

    data = (double **)malloc(NCM(n,m) * sizeof(double *));
    for (i=0; i<n; i++)
    {
        data[i] = (double *)malloc(n * sizeof(double));
    }

    printf("Enter Constraint equations with < , > or =\n");

    A = (double **)malloc(m * sizeof(double *));
    for (i=0; i<m; i++)
    {
        A[i] = (double *)malloc(n * sizeof(double));
    }
    A_BFS = (double **)malloc(m * sizeof(double *));
    for (i=0; i<m; i++)
    {
        A_BFS[i] = (double *)malloc((n+m) * sizeof(double));
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
            A_BFS[i][j] = A[i][j]; 
    	}

        scanf("%c",&dummy);
        scanf("%c",&symbol);

        if(symbol == '<')
        {
            for(j=n;j<n+m;j++)
            {
                A_BFS[i][j] =0;
            }
            A_BFS[i][i+n] =1;
            additional[0][i] = 1;
            additional[1][i] = 0;
            totalVariables++;
        }
        else if(symbol == '>')
        {
            for(j=n;j<n+m;j++)
            {
                A_BFS[i][j] =0;
            }
            A_BFS[i][i+n] = -1;
            additional[0][i] = -1;
            additional[1][i] = 1;
            totalVariables++;
            totalVariables++;
        }
        else
        {
            for(j=n;j<n+m;j++)
            {
                A_BFS[i][j] =0;
            }
            A_BFS[i][i+n] = 0;
            additional[0][i] = 0;
            additional[1][i] = 1;
            totalVariables++;
        }
        scanf("%lf",&B[i]);
    }

    tableau = (double **)malloc((m+2) * sizeof(double *));
    for (i=0; i<m+2; i++)
    {
        tableau[i] = (double *)malloc((totalVariables+3) * sizeof(double));
    }
    dummyTableau = (double **)malloc((m+2) * sizeof(double *));
    for (i=0; i<m+2; i++)
    {
        dummyTableau[i] = (double *)malloc((totalVariables+3) * sizeof(double));
    }

    
    combination = (int *)malloc(m*sizeof(int));
    expression = (double *)malloc((n+1)*sizeof(double));

    printf("Enter the coefficients of expression and 1 if you want to maximise else -1\n");
    for(i=0;i<=n;i++)
    {
        scanf("%lf",&expression[i]);
    }

    generateSolution(combination,A_BFS,B,0,0,m,n+m);
    makeTableau(A,B,expression,m+2,totalVariables+3,tableau,additional,n);
    totalIterations = applySimplex(tableau,m+2,totalVariables+3);

    printf("----------------------------\n");
    printf("Select among the following options\n");
    printf("1 List all BFS\n");
    printf("2 Number of iterations to solve the problem\n");
    printf("3 List of all Non-basic variables along with net evaluations in ith iteration\n");
    printf("4 List of Basic variables along with min ratios in ith iteration\n");
    printf("5 simplex table of ith iteration\n");
    printf("6 optimal solution (if exists else a report on others)\n");
    printf("0 to exit\n");
    while(1)
    {
        printf("----------------------------\n");
        printf("option: ");
        scanf("%d",&input);
        printf("----------------------------\n");
        if(input==1)
        {
            printf("Printing all BFS\n");
            for(i=0;i<count;i++)
            {
                flag = 0;
                for(j=0;j<n;j++)
                {
                    if(data[i][j]<0) flag=1;
                }
                if(flag==0)
                {
                    printf("[ ");
                    for(j=0;j<n;j++)
                    {
                        printf("%.3lf ", data[i][j]);
                    }
                    printf("]");
                    printf("\n");
                }
            }
        }
        else if(input==2)
        {
            if(unbounded==1)
            {
                printf("The expression is unbounded\n");
                printf("Total iterations performed before realising= %d\n",totalIterations);
                
            }
            else if(infeasible==1)
            {
                printf("The LPP is infeasible\n");
                printf("Total iterations performed before realising= %d\n",totalIterations);
            } 
            else
            {
                printf("Total iterations= %d\n",totalIterations);
            }
        }
        else if(input==3)
        {
            makeTableau(A,B,expression,m+2,totalVariables+3,dummyTableau,additional,n);
            printf("Enter the value of i: ");
            scanf("%d",&i);
            if(i<=totalIterations)
            {
                printf("Non Basic Variables and their net evaluations for %d-th iteration are:\n",i);
                while(i--)
                {
                    performSimplexIteration(dummyTableau,m+2,totalVariables+3,-1);
                }
                for(i=2;i<totalVariables+2;i++)
                {
                    for(j=2;j<m+2;j++)
                    {
                        if(dummyTableau[0][i]==dummyTableau[j][0]) break;
                    }
                    if(j==m+2)
                    {
                        temp = ((int)dummyTableau[0][i])%100;
                        if(dummyTableau[0][i] > 300) printf("artif%d ",temp); 
                        else if(dummyTableau[0][i] > 200) printf("slack%d ",temp);
                        else if(dummyTableau[0][i]>100) printf("surplus%d ",temp);
                        else printf("X%d ",temp);
                        printf("%.3lf\n",dummyTableau[1][i]);
                    }
                }
                printf("\n");
            }
            else printf("i cannot be more than the totalIterations (= %d )\n",totalIterations);
        }
        else if(input==4)
        {
            makeTableau(A,B,expression,m+2,totalVariables+3,dummyTableau,additional,n);
            printf("Enter the value of i: ");
            scanf("%d",&i);
            if(i<=totalIterations)
            {
                printf("Basic Variables along with ratios for %d-th iteration are: ",i);
                while(i--)
                {
                    performSimplexIteration(dummyTableau,m+2,totalVariables+3,-1);
                }
                double minNeg = 0;
                int colIdx;
                for(i=2;i<m+1;i++)
                {
                    if(dummyTableau[1][i]<minNeg)
                    {
                        minNeg = dummyTableau[1][i];
                        colIdx = i;
                    }
                }
                for(j=2;j<m+2;j++)
                {
                    temp = ((int)dummyTableau[j][0])%100;
                    if(dummyTableau[j][0] > 300) printf("artif%d ",temp); 
                    else if(dummyTableau[j][0] > 200) printf("slack%d ",temp);
                    else if(dummyTableau[j][0]>100) printf("surplus%d ",temp);
                    else printf("X%d ",temp); 

                    if(dummyTableau[j][colIdx]>0) printf(" %.3lf/%.3lf = %.3lf\n",dummyTableau[j][totalVariables+2],dummyTableau[j][colIdx],(dummyTableau[j][totalVariables+2]/dummyTableau[j][colIdx]));
                    else printf("non positive entering column variable");
                }
                printf("\n"); 
            }
            else printf("i cannot be more than the totalIterations (= %d )\n",totalIterations);
        }
        else if(input==5)
        {
           makeTableau(A,B,expression,m+2,totalVariables+3,dummyTableau,additional,n);
            printf("Enter the value of i: ");
            scanf("%d",&i);
            if(i<=totalIterations)
            {
                printf("Simplex table of %d-th iteration\n",i);
                while(i--)
                {
                    performSimplexIteration(dummyTableau,m+2,totalVariables+3,-1);
                }
                printTableau(dummyTableau,m+2,totalVariables+3);
            } 
        }
        else if(input==6)
        {
            if(unbounded==1) printf("The expression is unbounded\n");
            else if(infeasible==1) printf("The LPP is infeasible\n");
            else 
            {
                if(expression[n]==1) printf("The maximum value of the expression ");
                else printf("The minimum value of the expression ");
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
                        if(i == (int)tableau[j][0])
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
            }
        }
        else
        {
            break;
        }
    }
	return 0;
}
