#include<stdio.h>
#include<stdlib.h>

double** MatrixMultiplication(double **M,double **N,int n)
{
  int i,j,k;
  double **R;
  R=(double **)malloc(n*sizeof(double *));
  for(i=0;i<n;i++) 
    R[i]=(double*)malloc(n*sizeof(double));
  for(i=0;i<n;i++) 
    for(j=0;j<n;j++)
      {
	R[i][j]=0;
	for(k=0;k<n;k++)
	  R[i][j]+=M[i][k]*N[k][j];
      }
  return R;
}

void VectorMultiply(double *C,double **M,double *C1,int n)
{
  int i,j;
  for(i=0;i<n;i++)
    {
      C[i]=0;
      for(j=0;j<n;j++)
	C[i]+=M[i][j]*C1[j];
    }
}


void inverse(double **A,double **B,int n)
{
  double **M,*C,**E,*C1,temp;
  int i,j;
  M=(double **)malloc(n*sizeof(double *));
  for(i=0;i<n;i++) 
    M[i]=(double*)malloc(n*sizeof(double));
  E=(double **)malloc(n*sizeof(double *));
  for(i=0;i<n;i++) 
    E[i]=(double*)malloc(n*sizeof(double));
  C=(double*)malloc(n*sizeof(double));
  C1=(double*)malloc(n*sizeof(double));
  for(i=0;i<n;i++) 			//E=M=I
    for(j=0;j<n;j++)
      if(i==j)
	{
	  E[i][j]=1;
	  M[i][j]=1;
	}
      else
	{
	  E[i][j]=0;
	  M[i][j]=0;
	}
  for(i=0;i<n;i++)
    { 			//Iterations on columns
      for(j=0;j<n;j++)
	C1[j]=A[j][i];	//ith column
      VectorMultiply(C,M,C1,n);		//Multiply M*ith column
      temp=C[i]; 				
      for(j=0;j<n;j++)
	{		//Compute neta
	  if(j==i) C[j]=1/temp;	
	  else C[j]=(-1)*(C[j]/temp);
	}
      for(j=0;j<n;j++)
	{			//Ei: make (i-1)th column = I , ith row=neta
	  if(i!=0)
	    {
	      if((i-1)==j)
		E[j][i-1]=1;
	      else
		E[j][i-1]=0;
	    } 
	  E[j][i]=C[j];
	}
      M=MatrixMultiplication(E,M,n);	
    }
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      B[i][j]=M[i][j];		//M after all iterations is the resulting inverse.
}


void ShowMat(double** A,int m,int n)
{
  int i,j;
  for(i=0;i<m;i++)
    {
      for(j=0;j<n;j++)
	printf("%lf\t",A[i][j]);			
      printf("\n");
    }
}

void Mult(double **R,double **A,double **B,int a,int b,int c)
{
  int i,j,k;
  for(i=0;i<a;i++)
    {
      for(j=0;j<c;j++)
	{
	  R[i][j]=0;
	  for(k=0;k<b;k++)
	    R[i][j]+=A[i][k]*B[k][j];
	}
    }
}

int Entering(double **A,double **M,double*Z,int *T,int m,int n)
{
  double **Y;
  int i,j;
  Y=(double **)malloc(sizeof(double*));
  Y[0]=(double *)malloc(m*sizeof(double));
  double **Cb;
  Cb=(double**)malloc(sizeof(double*));
  Cb[0]=(double*)malloc(m*sizeof(double));
  for(i=0;i<m;i++)
    Cb[0][i]=Z[T[i+n]];
  Mult(Y,Cb,M,1,m,m);
  double min=1;
  int index=-1;
  for(i=0;i<n;i++)
    {
      double sum=0;
      for(j=0;j<m;j++)
	sum+=A[j][T[i]]*Y[0][j];
      sum-=Z[T[i]];
      if((sum<0)&&(sum<min))
	{
	  sum=min;
	  index=i;
	}
    }
  return index;
}

int Departing(double **A,double **B,double **M1,double **Xb,int *T,int m,int n,int e)
{
  int i;
  Mult(Xb,M1,B,m,m,1);
  double **Alj,**Pj;
  Alj=(double**)malloc(m*sizeof(double*));
  for(i=0;i<m;i++)
    Alj[i]=(double*)malloc(sizeof(double));
  Pj=(double**)malloc(m*sizeof(double*));
  for(i=0;i<m;i++)
    {
      Pj[i]=(double*)malloc(sizeof(double));
      Pj[i][0]=A[i][T[e]];
    }
  Mult(Alj,M1,Pj,m,m,1);
  int index=-1;
  double min=0,sum;
  for(i=0;i<m;i++)
    if((Alj[i][0]<=0)||(Xb[i][0]<=0))
      continue;
    else
      {
	sum=Xb[i][0]/Alj[i][0];
	if((index==-1)||(sum<min))
	  {
	    min=sum;
	    index=i;
	  }
      }
  return index;
}

void Result(double *Z,int *T,double **Xb,int m,int n)
{
  double r=0;
  int i, f, j;
  for(i=0;i<m;i++)
    r+=Xb[i][0]*Z[T[n+i]];
  printf("\nOptimum solution is Z = %lf.\n",r);
  printf("It occurs for :\n");
  for(i=0;i<m;i++)
    {
      f=0;
      printf("x[%d] = ", i+1);
      for(j=0; j<m; j++)
	if(T[n+j]==i)
	  {
	    printf("%lf\n", Xb[i][0]);
	    f=1;
	  }
      if(f!=1)
	printf("0\n");
    }
}

void Revised(double**A,double**B,double*Z,int m,int n)
{
  int *T;						//T= Matrix storing positions of NBV and BV
  int i,j;
  T=(int*)malloc((n+m)*sizeof(int));
  for(i=0;i<n+m;i++)
    T[i]=i;
  double **M,**M1;						//M=Basis Matrix B
  M=(double**)malloc(m*sizeof(double*));
  for(i=0;i<m;i++)
    M[i]=(double*)malloc(m*sizeof(double));
  M1=(double**)malloc(m*sizeof(double*));
  for(i=0;i<m;i++)
    M1[i]=(double*)malloc(m*sizeof(double));
  for(i=0;i<m;i++)
    for(j=0;j<m;j++)
      if(i==j)
	M[i][j]=1;
      else
	M[i][j]=0;
  double **Xb;
  Xb=(double**)malloc(m*sizeof(double*));
  for(i=0;i<m;i++)
    Xb[i]=(double*)malloc(sizeof(double));
  while(1)
    {
      int e;
      //printf("\n");
      //ShowMat(M,m,m);
      inverse(M,M1,m);
      //printf("\n");
      //ShowMat(M1,m,m);
      e=Entering(A,M1,Z,T,m,n);			//Gives departing variable
      if (e<0)
	{
	  Mult(Xb,M1,B,m,m,1);
	  Result(Z,T,Xb,m,n);
	  return;
	}
      int d;
      d=Departing(A,B,M1,Xb,T,m,n,e);
      if(d<0)
	{
	  printf("\nInfeasible\n");
	  return;
	}
      for(i=0;i<m;i++)
	M[i][d]=A[i][T[e]];
      d+=n;
      int swap;
      swap=T[d];
      T[d]=T[e];
      T[e]=swap;
    }
}

int main()
{
  int n,m,i,j;
  printf("Max :       \tZ = cx\nSubject to :\tAx<=B\n\t\tx>=0\n");
  printf("Enter the number of variables : ");
  scanf("%d",&n);
  printf("\nEnter the number of constraints(ignore positivity of variables) : ");
  scanf("%d",&m);
  printf("\nEnter the coefficients of variables in the objective function :\n");
  double *Z;
  Z=(double*)malloc((n+m)*sizeof(double));
  for(i=0;i<n;i++)
    {
      printf("x[%d] : ", i+1);
      scanf("%lf",&Z[i]);
    }
  for(i=n;i<n+m;i++)
    Z[i]=0;
  printf("\na(i,1)x1 + a(i,2)x2 + ... +a(i,n)xn <= bi\n");
  printf("Enter the coefficients of variables in constraints :\n");
  double **A;
  A=(double**)malloc(m*sizeof(double*));
  double **B;
  B=(double**)malloc(m*sizeof(double*));
  for(i=0;i<m;i++)
    B[i]=(double *)malloc(sizeof(double));
  for(i=0;i<m;i++)
    A[i]=(double *)malloc((n+m)*sizeof(double));
  for(i=0;i<m;i++)
    {
      for(j=0;j<=n;j++)
	if(j!=n)
	  {
	    printf("a(%d, %d) = ", i+1, j+1);
	    scanf("%lf",&A[i][j]);
	  }
	else
	  {
	    printf("b%d = ", i+1);
	    scanf("%lf",&B[i][0]);
	    printf("\n");
	  }
      for(j=n;j<n+m;j++)
	if(i==(j-n))
	  A[i][j]=1;
	else
	  A[i][j]=0;
    }
  Revised(A,B,Z,m,n);
}
