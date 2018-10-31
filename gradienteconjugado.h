
/**
 * @file gradienteconjugado.h
 * @author GRR20176143 Cláudio Torres Júnior
 * @author GRR20171607 Gabriela Stein
 * @date 16 Sep 2018
 * @brief Arquivo de header para as funções e definições de tipo.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h> /* memset */

// Parâmetros para teste de convergência

/**
* @brief Valor de epsilon para teste de convergência
*/

#define EPS 1.0e-4

/**
* @brief Estrutura que contêm todos os parâmetros passados na execução do programa
*/

typedef struct parametro{
	long int n, //dimensao
			 k, //n de diagonais
		     i; //max de iteracao
	double p, //pré-condicionador
		   e; //erro aproximado
	char *o;
	int op; //nos diz se o erro eh opcional
}parametro;

/**
* @brief Estrutura que contêm o tempo inicial, final e a diferença entre os diagonais
*/

typedef struct tempo{
	double ini;
	double fim;
	double dif;
}tempo;

void multMatMat(double *pri, double *sec, long int dgn, long int tam, double *mult);
void multMatVet(double *pri, double *sec, long int inicio, long int dgn, long int tam, double *mult);
double multVetVet(double *pri, double *sec, long int dgn, long int tam);
void trasformaSistema(double *A, double *B, double *Atf, double *Btf, parametro par);
void transposta(double *A, double *T, parametro par);
void preCondicionador(double *M, double *A, parametro par);
int gradienteConjugado(double *A, double *B, parametro par);
void criaMatrizes(double *A, double *L, double *U, double *D, long int tam);
double maxVetor(double *V, parametro par);
void Cholesky(double *M, long int tam);
void liberaVet(double *M, double *X, double *r, double *v, double *z, 
	double *y, double *Xant, double *erroAproximadoA, double *erroIt);