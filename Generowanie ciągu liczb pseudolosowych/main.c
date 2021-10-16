#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double gen_1(int a, int c, int m_p)
{
  static long int x = 10;
  long int m = pow(2, m_p);

  x = (a * x + c) % m;

  return (double)x / ((double)m + 1.0);
}

double F(double x, double mi, double d)
{
  if (x <= mi)
  {
      return (-(1.0 / (d * d))) * (-(x * x) / 2.0 + mi * x) + x / d;
  }
  return (-(1.0 / (d * d))) * ((x * x) / 2.0 - mi * x + mi * mi) + x / d;
}

int main()
{
  FILE *file_U = fopen("U.dat", "w");
  FILE *file_jednorodny = fopen("jednorodny.dat", "w");
  FILE *file_trojkatny = fopen("trojkatny.dat", "w");
 
  int N = 10000, K = 12, j;
  double x1[N], x2[N];
  double sum_avr = 0., sum_stdev = 0.;
  double avr, stdev, dx = 1.0 / K;
  int nj1[K], nj2[K];

  for (int i = 0; i < K; i++)
  {
    nj1[i] = 0;
    nj2[i] = 0;
  }
  int a = 123, c = 1, m_p = 15;

  for (int i = 0; i < N; i++)
  {
    x1[i] = gen_1(a, c, m_p);
    if (i > 0)
    {
      fprintf(file_U, "%f\t%f\n", x1[i - 1], x1[i]);
    }
    sum_avr += x1[i - 1];
    j = (int)(x1[i] / dx);
    nj1[j]++;
  }
  fprintf(file_U, "\n\n");

  avr = sum_avr / N;

  for (int i = 0; i < N; i++)
  {
    sum_stdev += (x1[i] - avr) * (x1[i] - avr);
  }
  stdev = sqrt(sum_stdev / N);

  for (int j = 0; j < K; j++)
  {
    fprintf(file_jednorodny, "%f\t%f\n", (j + 0.5) * dx, (double)nj1[j] / N);
  }
  fprintf(file_jednorodny, "\n\n");

  a = 69069, c = 1, m_p = 32;
  sum_avr = 0., sum_stdev = 0.;

  for (int i = 0; i < N; i++)
  {
    x2[i] = gen_1(a, c, m_p);
    if (i > 0)
    {
      fprintf(file_U, "%f\t%f\n", x2[i - 1], x2[i]);
    }
    sum_avr += x2[i - 1];
    j = (int)(x2[i] / dx);
    nj2[j]++;
  }
  fprintf(file_U, "\n\n");

  avr = sum_avr / N;

  for (int i = 0; i < N; i++)
  {
    sum_stdev += (x2[i] - avr) * (x2[i] - avr);
  }
  stdev = sqrt(sum_stdev / N);

  for (int j = 0; j < K; j++)
  {
    fprintf(file_jednorodny, "%f\t%f\n", (j + 0.5) * dx, (double)nj2[j] / N);
  }

  N = 1000;
  K = 10;
  double x3[N], pj[K];
  int nj3[K];
  double mi = 4.0, d = 3.0, ksi_1, ksi_2;
  dx = ((mi + d) - (mi - d)) / K;

  for (int i = 0; i < K; i++)
  {
    nj3[i] = 0;
  }

  for (int i = 0; i < N; i++)
  {
    ksi_1 = gen_1(a, c, m_p); 
    ksi_2 = gen_1(a, c, m_p);
    x3[i] = mi + (ksi_1 + ksi_2 - 1.0) * d;

    j = (int)((x3[i] - (mi - d)) / dx);
    nj3[j]++;
  }

  for (int j = 0; j < K; j++)
  {
    pj[j] = F(((mi - d) + (j + 1) * dx), mi, d) - F(((mi - d) + j * dx), mi, d);

    fprintf(file_trojkatny, "%f\t%f\t%f\n", (mi - d) + (j + 0.5) * dx, (double)nj3[j] / N, pj[j]);
  }
  double chi2 = 0.0; 

  for (int j = 0; j < K; j++)
  {
    chi2 += pow((nj3[j] - N * pj[j]), 2.0) / (N * pj[j]);
  }
  fclose(file_U);
  fclose(file_jednorodny);
  fclose(file_trojkatny);

  return 0;
}