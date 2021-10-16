#include <stdio.h>
#include <math.h>

#define N 7 
#define IT_MAX 12

void matrixVector(double (*A)[N], double * v, double * result, int size) {
	double temp;
		for (int i = 0; i < size; i++) {
			temp = 0.0;
			for (int j = 0; j < size; j++) 
				temp += A[i][j] * v[j];
			result[i] = temp;
		}
}

double vectorVector(double * v1, double * v2, int size) {
	double result = 0.0;

	for (int i = 0; i < size; i++)
		result += v1[i] * v2[i];

	return result;
}

void matrixMatrix(double (*A)[N], double (*B)[N], double (* result)[N], int size) {
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			result[i][j] = 0.0;
			for (int k = 0; k < size; k++)
				result[i][j] += A[i][k] * B[k][j];
		}
	}
}

void orthogonalization(double * v1, double * v2, int size) {
	double temp = vectorVector(v1, v2, size);

	for (int i = 0; i < size; i++)
		v1[i] = v1[i] - temp * v2[i];
}

int main() 
{
	double A[N][N], X[N][N], D[N][N], przyblizenia[IT_MAX][N];
	double x[N], next[N], lambda = 0.0, norma;
	int K_ww = N;
	FILE * f1 = fopen("lambda.dat", "w");
	FILE * f2 = fopen("macierz_X.dat", "w");
	FILE * f3 = fopen("macierz_D.dat", "w");
	
	//printf("Macierz A:\n");
	for(int i = 0; i < N; i ++) 
	{
		for(int j = 0; j < N; j ++) 
		{
			A[i][j] = 1.0 / sqrt(2.0 + fabs(i-j));
			//printf("%g ", A[i][j]);
		}
		//printf("\n");
	}
	//printf("\n");

	for(int k = 0; k < K_ww; k++) 
	{
		for(int i = 0; i < N; i++)
		{
			x[i] = 1.0;
		}

		for(int i = 0; i < IT_MAX; i++) 
		{
			matrixVector(A, x, next, N);

			for(int j = 0; j < k; j++) 
			{
				orthogonalization(next, X[j], N);
			}
			lambda =  vectorVector(next, x, N) / vectorVector(x, x, N);
			przyblizenia[i][k] = lambda;
			norma = sqrt(vectorVector(next, next, N));

			for(int j = 0; j < N; j++) 
			{
				x[j] = next[j]/norma;
			}
		}
		fprintf(f1, "%g\n", lambda);

		for (int i = 0; i < N; i++)
		{
			X[k][i] = x[i];
		}
	}

	for (int i = 0; i < N; i++) 
	{
		for (int j = 0; j < N; j++) 
		{
			fprintf(f2, "%12g ", X[j][i]);
		}
		fprintf(f2, "\n");
	}

	double tempMat[N][N];
	double tempVal;
	matrixMatrix(X, A, tempMat, N);

	for (int i = 0; i < N; i++) 
	{
		for (int j = i + 1; j < N; j++) 
		{
			tempVal = X[i][j];
			X[i][j] = X[j][i];
			X[j][i] = tempVal;
		} 
	}
	matrixMatrix(tempMat, X, D, N);

	for (int i = 0; i < N; i++) 
	{
		for (int j = 0; j < N; j++) 
		{
			fprintf(f3, "%15g ", D[j][i]);
		}
		fprintf(f3, "\n");
	}
	fclose(f1);
	fclose(f2);
	fclose(f3);

	return 0;
}