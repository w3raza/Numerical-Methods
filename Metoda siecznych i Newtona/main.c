#include<stdio.h>
#include<math.h>

double g1(double x)           
{ 
    return sin(x);        
}

double g2(double x)           
{ 
    return (x*x)/8;      
}

double f(double x)            
{ 
    return g1(x) - g2(x); 
}

double derivative_f(double x) 
{ 
    return cos(x)-(x/4);  
}


void save_function(const char* filename, double (*fun)(double), double start, double end)
{
    FILE* fp = fopen(filename, "w");

    for(double i = start; i<end; i+=0.1)
    {
        fprintf(fp, "%10g %10g\n", i, fun(i));
    }
    fclose(fp);
}

void newton(const char* filename, int it_max, double x0)
{
    FILE* fp = fopen(filename, "w");
    
    for(int i=1; i<=it_max; i++)
    {
        x0 = x0 - (f(x0)/derivative_f(x0));
        fprintf(fp, "%d %g %10g %10g\n", i, x0, f(x0), derivative_f(x0));
    }
}

void sieczne(const char* filename, int it_max, double a, double b)
{
    FILE* fp = fopen(filename, "w");

    for(int i=1; i<=it_max; i++)
    {
        double x = b - (f(b)*(b-a))/(f(b)-f(a));
        a = b;
        b = x; 
        fprintf(fp, "%d %g %g %g\n", i, x, f(b), f(a));
    }
    fclose(fp);
}

int main()
{
    double x = -8.0;
    double x0 = 8.0;
    double x1 = -8.1;
    double x2 = 8.1;

    save_function("g1.dat", g1, x, x0);
    save_function("g2.dat", g2, x, x0);
    save_function("f.dat", f, x, x0);

    newton("newton_ujemne.dat", 10, x);
    newton("newton_dodatnie.dat", 10, x0);

    sieczne("sieczne_ujemne.dat", 15, x, x1);
    sieczne("sieczne_dodatnie.dat", 15, x0, x2);

    return 0;
}
