#include <gsl/gsl_linalg.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define x_max 5
#define x_min -5

//f1(x)=1/(1+x^2)
double f_1(double x)
{
    return 1.0/(1.0 + x*x);
}
//f2(x)=cos(2x)
double f_2(double x)
{
    return cos(2.0*x);
}
double deriv(double x, double dx)
{
    return (f_1(x - dx) - 2*f_1(x) + f_1(x + dx)) / pow(dx, 2);
}

void wyznacz_M(double *xw, double *yx, gsl_vector *m, int n, double alfa, double beta);
double wyznacz_Sx(double *xw, double *yx, gsl_vector *m, int n, double x);

void printVector(gsl_vector *v, int n);
void printVectorAcc(gsl_vector *v, int n);
void printMatrix(gsl_matrix *M, int n);

int main()
{
    
    FILE *f1;
    f1 = fopen("f1.dat","w");
    if(f1 == NULL)
        exit(1);
    FILE *fdrv;

    FILE *f2;    
    f2 = fopen("f2.dat","w");
    if(f2 == NULL)
       exit(1);
 
    fdrv = fopen("pochodne.dat","w");
    if(fdrv == NULL)
    {
        printf("Error!");   
        exit(1);             
    }
    double alfa = 0;
    double beta = 0;
    int n = 5;
    double xw[n];
    double yx[n];
    double delta = ((double)x_max - (double)x_min) / ((double)n - 1.0);
    int i,f;
    double dx = 0.01;

    for(i=0;i<n;i++)
        xw[i] = (double)x_min + delta*i;

    for(i=0;i<n;i++)
    {
        yx[i] = f_2(xw[i]);
    }
    gsl_vector *m = gsl_vector_calloc(n);
    
    wyznacz_M(xw, yx, m, n, alfa, beta);

    for(double x = (double)x_min; x<= (double)x_max; x+= 0.1) 
    { 
        double temp = wyznacz_Sx(xw, yx, m, n, x);
        fprintf(f1,"%f %f\n",x, temp);
    }
    fprintf(f1, "\n\n");
    
    if(f == 1 && n == 5)
    {   
        double temp;
        for(i=0;i<n;i++)
        {
            temp = deriv(xw[i] , dx);
            fprintf(fdrv,"%f %f %f\n", xw[i], gsl_vector_get(m, i), temp);
        }
    } 

    gsl_vector_free(m);
    fclose(fdrv);
    
    return 0;
}

void wyznacz_M(double *xw, double *yx, gsl_vector *m, int n, double alfa, double beta)
{
    gsl_matrix * A = gsl_matrix_alloc(n, n);
    gsl_vector * d = gsl_vector_alloc(n);

    double h_i = ((double)x_max - (double)x_min) / ((double)n - 1.0);
    double lambda = h_i / (h_i + h_i);
    double mu = 1.0 - lambda;
    
    gsl_vector_set(d, 0, alfa);
    gsl_vector_set(d, n-1, beta);

    double value;
    for(int i = 1; i < n-1; i++)
    {
        value = (6.0 / (h_i + h_i)) * ( (yx[i+1] - yx[i]) / h_i - (yx[i] - yx[i-1]) / h_i);
        gsl_vector_set(d, i, value);
    }
    
    gsl_matrix_set_zero(A);

    for(int i = 0; i < n; i++)
    {
        gsl_matrix_set(A, i, i, 2);  
    }
    gsl_matrix_set(A, 0, 0, 1);
    gsl_matrix_set(A, n-1, n-1, 1);

   
    for(int i = 1; i < n-1; i++)
    {
        for(int j = 1; j < n-1 ; j++)
        {
            if(i == j)
            {
                gsl_matrix_set(A, i, j-1, mu);
                gsl_matrix_set(A, i, j+1, lambda);
            }
        }   
    }
    gsl_linalg_HH_svx( A, d);
    for(int i=0; i<n;i++)
    {
        gsl_vector_set(m, i, gsl_vector_get(d, i));
    }
    gsl_matrix_free(A); 
    gsl_vector_free(d);

}

double wyznacz_Sx(double *xw, double *yx, gsl_vector *m, int n, double x)
{
    int przedzial;
    double Sx;
    double h_i = ((double)x_max - (double)x_min) / ((double)n - 1.0);

    for(int i=1; i<n; i++)
    {
        if(xw[i-1] <= x && x <= xw[i])
        {
            przedzial = i-1;
            break;
        }
    }
    double A_i;
    double B_i;
    
    double tmp1 = gsl_vector_get(m , przedzial + 1);
    double tmp2 = gsl_vector_get(m , przedzial);
    A_i =( (yx[przedzial+1] - yx[przedzial]) / h_i ) - (h_i/6.0) * ( tmp1 - tmp2);    

    tmp1 = gsl_vector_get(m , przedzial);
    B_i = yx[przedzial] - tmp1 * ( (pow(h_i,2)) / 6.0 );
    Sx = gsl_vector_get(m , przedzial);
    Sx *= (pow((xw[przedzial+1] - x), 3) / (6.0*h_i));
    Sx += gsl_vector_get(m , przedzial + 1) * (pow((x - xw[przedzial]), 3) / (6.0*h_i));
    Sx += A_i * (x - xw[przedzial]);
    Sx += B_i;

    return Sx;
}

void printVector(gsl_vector *v, int n)
{
    puts("");
    double temp;
        for(int i=0; i<n ; i++)
        {
            temp =  gsl_vector_get(v , i);
            printf("| %.7f |\n", temp);     
        }
}

void printVectorAcc(gsl_vector *v, int n)
{
    puts("");
    double temp;
        for(int i=0; i<n ; i++)
        {
            temp =  gsl_vector_get(v , i);
            printf("| %.7e |\n", temp);     
        }
}

void printMatrix(gsl_matrix *M, int n)
{
    double temp;
    for(int i=0; i<n ; i++)
        {
            puts("");
            for(int j=0;j<n;j++)
            {
                temp =  gsl_matrix_get(M , i, j);
                printf("Value for [%d][%d] = %2.2f\n",i,j, temp); 
            }
        }
    puts("");
}