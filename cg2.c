#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "sistemarandom.h"

int main (int argc, char *argv[])
{
	double *A, *B;
	int i, j;
	srand(20182);
 	double num = 1.0;

	int k = atoi(argv[1]);
	int n = atoi(argv[2]);

	A = (double*)malloc(n*n*sizeof(double));
	B = (double*)malloc(n*sizeof(double));

	for(i = 0; i < n; i++){
		for(j = 0; j < n; j++){
			if(abs(i - j) > k/2){
				A[i*n + j] = 0;
			}else{
				A[i*n + j] = num++;//generateRandomA(i, j, k);
			}
		}
	}

	for(i = 0; i < n; i++){
		B[i] = generateRandomB(k);
	}

	/*---------------------------------------------------------------*/
	int col;
	int numZeros;
	int tam = n*k;
	int c;
	num = 1.0;

	double *vet = (double*)malloc(tam*sizeof(double));
	double *t = (double*)malloc(tam*sizeof(double));

	memset(vet, 0, sizeof(vet));
	memset(t, 0, sizeof(t));

	numZeros = k/2;
	col = k - k/2;
	c = 0;
	for(i = 0; i < n; i++){
		for(j = c; (j < col) && (j < n); j++){
			vet[i*k + (k/2 - i) + j] = num++;
		}
		col++;
		numZeros--;
		if(numZeros < 0){
			c++;
		}
	}

	numZeros = k/2;
	col = k - k/2;
	c = 0;
	for(i = 0; i < n; i++){
		for(j = c; (j < col) && (j < n); j++){
			t[i*k + (k/2 - i) + j] = vet[j*k + (k/2 - j) + i];
		}
		col++;
		numZeros--;
		if(numZeros < 0){
			c++;
		}
	}
	return(0);
}