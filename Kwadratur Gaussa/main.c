#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include"nrutil.h"
#include"nrutil.c"
#include"gauher.c"
#include"gauleg.c"
#include"gaulag.c"
#include"gammln.c"

double f1(double x)  { return log(x); }
double f2a(double x) { return pow(x-10, 2)*sin(4*x); }
double f2b(double x) { return f2a(x)*exp(-x); }
double f3a(double x) { return pow(x, 7)*pow(2, -x*x+x+4); }
double f3b(double x) { return f3a(x)*exp(-x*x); }

void zad1();
void zad2();
void zad3();

void pf1();
void pf2();
void pf3();

int main(void)
{
	zad1();
	zad2();
	zad3();

	pf1();
	pf2();
	pf3();

	return 0;
}

void zad1()
{
	FILE *fp = fopen("zad1.dat", "w");
	double teo = 10*log(10)-10;

	for(int n=5; n<=70; n++)
	{
		float* x = calloc(n+1, sizeof(float));
		float* w = calloc(n+1, sizeof(float));

		gauleg(0, 10, x, w, n);

		double sum = 0;

		for(int i=1; i<=n; i++)
			sum += w[i]*f1(x[i]);

		fprintf(fp, "%d %g\n", n, fabs(teo - sum));

		free(x);
		free(w);
	}

	fclose(fp);
}

void zad2()
{
	FILE *fp = fopen("zad2.dat", "w");

	for(int n=5; n<=70; n++)
	{
		float* x1 = calloc(n+1, sizeof(float));
		float* w1 = calloc(n+1, sizeof(float));

		float* x2 = calloc(n+1, sizeof(float));
		float* w2 = calloc(n+1, sizeof(float));

		gaulag(x1, w1, n, 0);
		gauleg(0, 10, x2, w2, n);

		double sum1 = 0;
		double sum2 = 0;

		for(int i=1; i<=n; i++)
		{
			sum1 += w1[i]*f2a(x1[i]);
			sum2 += w2[i]*f2b(x2[i]);
		}

		fprintf(fp, "%d %g %g\n", n, fabs(22.95461022 - sum1), fabs(22.95461022 - sum2));

		free(x1);
		free(w1);
		free(x2);
		free(w2);
	}

	fclose(fp);
}

void zad3()
{
	FILE *fp = fopen("zad3.dat", "w");

	for(int n=5; n<=70; n++)
	{
		float* x1 = calloc(n+1, sizeof(float));
		float* w1 = calloc(n+1, sizeof(float));

		float* x2 = calloc(n+1, sizeof(float));
		float* w2 = calloc(n+1, sizeof(float));

		gauher(x1, w1, n);
		gauleg(-10, 15, x2, w2, n);

		double sum1 = 0;
		double sum2 = 0;

		for(int i=1; i<=n; i++)
		{
			sum1 += w1[i]*f3a(x1[i]);
			sum2 += w2[i]*f3b(x2[i]);
		}

		fprintf(fp, "%d %g %g\n", n, fabs(14.83995751 - sum1), fabs(14.83995751 - sum2));

		free(x1);
		free(w1);
		free(x2);
		free(w2);
	}

	fclose(fp);
}

void pf1()
{
	FILE* fp = fopen("f1.dat", "w");

	for(double x=0; x<10; x += 0.1)
		fprintf(fp, "%g %g\n", x, f1(x));

	fclose(fp);
}

void pf2()
{
	FILE* fp = fopen("f2.dat", "w");

	for(double x=0; x<20; x+=0.1)
		fprintf(fp, "%g %g\n", x, f2b(x));

	fclose(fp);
}

void pf3()
{
	FILE* fp = fopen("f3.dat", "w");

	for(double x=-10; x<15; x+=0.1)
		fprintf(fp, "%g %g\n", x, f3b(x));

	fclose(fp);
}