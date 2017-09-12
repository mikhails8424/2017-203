#include "seskinmvs.h"

/**
 * Метод Гаусса
 */
void seskinmvs::lab1()
{
double t;
for (int k = 0; k < N; k++)
{
t = A[k][k];
for (int j = 0; j < N; j++)
A[k][j] = A[k][j] / t;
b[k] =b[k]/t;
for (int i = k + 1; i < N; i++)
{
t = A[i][k];
for (int j = 0; j < N; j++)
{
A[i][j] =A[i][j]- A[k][j] * t;
}
b[i] =b[i]- b[k] * t;
}
}
for (int k = N - 1; k > 0; k--)
{
for (int i = k - 1; i >= 0; i--)
{
t = A[i][k];
for (int j = 0; j < N; j++)
A[i][j] =A[i][j]- A[k][j] * t;
b[i] =b[i] - b[k] * t;
}
}
for(int i=0; i<N; i++)
x[i]=b[i];
}


/**
 * Метод Гаусса с выбором главного элемента
 */
void seskinmvs::lab2()
{
double t;
for (int k = 0; k < N; k++)
{
int maxelem=k;
for(int i=k+1;i<N;i++)
if(abs(A[i][k]) > abs(A[maxelem][k]))
maxelem=i;
for(int i=0;i<N;i++)
std::swap(A[k][i],A[maxelem][i]);
std::swap(b[k],b[maxelem]);
t = A[k][k];
for (int j = 0; j < N; j++)
A[k][j] = A[k][j] / t;
b[k] =b[k]/t;
for (int i = k + 1; i < N; i++)
{
t = A[i][k];
for (int j = 0; j < N; j++)
{
A[i][j] =A[i][j]- A[k][j] * t;
}
b[i] =b[i]- b[k] * t;
}
}
for (int k = N - 1; k > 0; k--)
{
for (int i = k - 1; i >= 0; i--)
{
t = A[i][k];
for (int j = 0; j < N; j++)
A[i][j] =A[i][j]- A[k][j] * t;
b[i] =b[i] - b[k] * t;
}
}

for(int i=0; i<N; i++)
x[i]=b[i];
}



/**
 * Метод квадратного корня (метод Холецкого)
 */
void seskinmvs::lab3()
{
double **L = new double*[N];
double **LT = new double*[N];
double *D = new double[N];
for(int i=0; i<N; i++)
{
L[i] = new double[N];
LT[i] = new double[N];
D[i]=0;
for(int j=0; j<N; j++)
{
L[i][j]=0;
}
}

for(int i=0; i<N; i++)
{
double isum = 0;
for(int k=0; k<i; k++)
isum += std::abs(L[k][i]) * std::abs(L[k][i]) * D[k];

if((A[i][i] - isum) > 0)
D[i] = 1;
else
D[i] = -1;

L[i][i] = sqrt(std::abs(A[i][i] - isum));

for(int j=i+1; j<N; j++)
{
double sum=0;
for(int k=0; k<i; k++)
sum += L[k][i] *
L[k][j] * D[k];
L[i][j] = (A[i][j] - sum) / (L[i][i] * D[i]);
}
}
for(int i=0; i<N; i++)
for(int j=0; j<N; j++)
LT[i][j] = L[j][i];
for(int i = 0; i < N; i++)
for(int j = 0; j < N; j++)
LT[i][j] *= D[j];

double *y = new double[N];
for(int i=0; i<N; i++)
{
double s=0;
for(int j=0; j<i; j++)
s += y[j] * LT[i][j];
y[i]=(b[i]-s)/LT[i][i];
}
//L * x = y
for(int i=N-1; i>=0; i--)
{
double s=0;
for(int j=i+1; j<N; j++)
s += x[j] * L[i][j];
x[i]=(y[i] - s)/L[i][i];
}

for(int i=0; i<N; i++)
{
delete [] L[i];
delete [] LT[i];
}
delete [] L;
delete [] LT;
delete [] D;
delete [] y;
}



/**
 * Метод прогонки
 */
void seskinmvs::lab4()
{

}



/**
 * Метод Якоби
 */
void seskinmvs::lab5()
{
double *xold = new double[N];
for (int i=0; i<N; i++)
{
x[i]=0; // первое решение
}
double p=0.0;
double eps=1e-20;
int k=0;
do
{
k++;
p=0.0;
for(int i=0; i<N; i++)
xold[i]=x[i]; //x old-старое решение
for(int i=0; i<N; i++)
{
double s=0;
for(int j=0; j<i; j++)
s += A[i][j] * xold[j];
for(int j=i+1; j<N; j++)
s += A[i][j] * xold[j];
x[i]=(b[i] - s)/A[i][i];
}
p= std::abs(xold[0]-x[0]);
for(int i=0; i<N; i++)
{
if(std::abs(xold[i]-x[i]) > p)
p = std::abs(xold[i]-x[i]);
}
} while(p >= eps);
std::cout « "Чиcло итераций : " « k « std::endl;

delete [] xold;
}



/**
 * Метод Зейделя
 */
void seskinmvs::lab6()
{

}



/**
 * Один из градиентных методов
 */
void seskinmvs::lab7()
{

}
