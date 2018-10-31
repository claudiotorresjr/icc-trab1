
/*! \mainpage Primeiro Trabalho de Introdução à Computação Científica 2018/2
 *
 * \section introducao Introdução
 *		O objetivo deste trabalho é implementar um programa computacional
 *		para testar a eficiência do método dos Gradientes Conjugados com
 *		diferentes pré-condicionadores para a resolução de sistemas esparsos
 *		nos quais a matriz A é k-diagonal.
 *
 * \subsection autores Autores
 * 		GRR20176143 Claúdio Torres Júnior  \n
 *		GRR20171607 Gabriela Stein
 *
 * \subsection professor Professor
 * 		Daniel Weingaertner
 * 
 */

/**
 * @file cgSolver.c
 * @author GRR20176143 Cláudio Torres Júnior
 * @author GRR20171607 Gabriela Stein
 * @date 16 Sep 2018
 * @brief Arquivo do programa principal.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
//#include "utils.h"
#include <unistd.h> // *POSIX* Para o getopt() original
#include <getopt.h> // *GNU* Para o getopt_long()

#include "sistemarandom.h"
#include "gradienteconjugado.h"

/**
 * @brief Função para explicar os parametros usados pelo programa, finaliza a execução.
*/

void ajuda(){
	printf("Argumentos passados incorretos. Deve ser no formato:\n");
	printf("cgSolver -n <n> -k <k> -p <ω> -i <i> -e <ε> -o <arquivo_saida>\n");
	printf("onde:");
	printf("<n> = obrigatorio. dimensao do sistema (n>10)\n");
	printf("<k> = obrigatorio. numero de diagonais. (k>1 e k impar)\n");
	printf("<ω> = obrigatorio. pre-condicionador\n");
	printf("\tω=0.0: sem pre-condicionador\n");
	printf("\t0.0 < ω < 1.0: pré-condicionador de Jacobi\n");
	printf("\tω=1.0 pré-condicionador de Gauss-Seidel\n");
	printf("\t1.0 < ω < 2.0: pré-condicionador SSOR\n");
	printf("<i> = obrigatorio. maximo de iteracoes\n");
	printf("<ε> = opcional. erro aproximado absoluto\n");
	printf("<arquivo_saida> = obrigatorio. caminho completo para o arquivo que vai conter a solução.\n");
	exit(-1);
}

/**
 * @brief Função para ler os parametros do programa
 * @param argc Contagem dos parâmetros
 * @param argv Vetor de strings contendo os parâmetros
 * @param par Struct usada para atribuição dos parâmetros
*/

void opcoes(int argc, char *argv[], parametro *par){
	int opt;

	while((opt = getopt(argc, argv, "n:k:p:i:e:o:")) != -1){
		switch(opt){
			case 'n':
				if(isdigit(*argv[optind - 1]))
					par->n = atoi(optarg);
				else{
					fprintf(stderr, "-n precisa de um parametro numerico positivo.\n");
					exit(-1);
				}
				if(par->n <= 10){
					fprintf(stderr, "n precisa ser > 10.\n");
					exit(-1);
				}
				break;
			case 'k':
				if(isdigit(*argv[optind - 1]))
					par->k = atoi(optarg);
				else{
					fprintf(stderr, "-k precisa de um valor associado numerico positivo.\n");
					exit(-1);
				}
				if( (par->k > par->n) || (par->k <= 1) || ((par->k % 2) == 0) ){
					fprintf(stderr, "k deve ser menor que n, ser > 1 e impar.\n");
					exit(-1);
				}
				break;
			case 'p':
				if(isdigit(*argv[optind - 1]))
					par->p = atof(optarg);
				else{
					fprintf(stderr, "-p precisa de um valor associado numerico positivo.\n");
					exit(-1);
				}
				if(par->p < 0.0 || par->p >= 2.0){
					fprintf(stderr, "p precisa estar: 0.0 <= p < 2.0\n");
					exit(-1);
				}
				break;
			case 'i':
				if(isdigit(*argv[optind - 1]))
					par->i = atoi(optarg);
				else{
					fprintf(stderr, "-i precisa de um valor associado numerico positivo.\n");
					exit(-1);
				}
				break;
			case 'e':
				if(strcmp(optarg, "-o") == 0){
					par->e = EPS;
					par->op = 1;
					par->o = argv[optind];
				}
				else{
					par->e = atof(optarg);
					par->op = 0;
					if (par->e <= 0.0){
						fprintf(stderr, "-e precisa de um valor positivo.\n");
						exit(-1);
					}
				}
				break;
			case 'o':
				par->o = optarg;
				break;
			case '?':
				fprintf(stderr, "----Parametro passado incorreto----\n");
				ajuda();
				break;
		}
	}  
}

/**
 * @brief Função principal do programa
 * @param argc Contagem dos parâmetros
 * @param argv Vetor de strings contendo os parâmetros
 * @return Retorna 0 caso tenha concluído com êxito.
*/

int main (int argc, char *argv[])
{
	long int i, j, col, numZeros, c;
	parametro par;

	srand(20182);

	if(argc != 13 && argc != 12){
		ajuda();
	}
	
	opcoes(argc, argv, &par);


	long int tam = par.n*par.k;
	numZeros = par.k/2;
	col = par.k - par.k/2;
	c = 0;


	double *A = (double*)malloc(tam*sizeof(double));
	double *B = (double*)malloc((par.n + numZeros)*sizeof(double));

	memset(A, 0, sizeof(*A));
	memset(B, 0, sizeof(*B));


	for(i = 0; i < par.n; ++i){
		for(j = c; (j < col) && (j < par.n); ++j){
			A[i*par.k + (par.k/2 - i) + j] = generateRandomA(i, j, par.k);
		}
		col++;
		numZeros--;
		if(numZeros < 0){
			c++;
		}
	}

	numZeros = par.k/2;
	for(i = numZeros; i < (par.n + numZeros); ++i){
		B[i] = generateRandomB(par.k);
	}

	//numZeros = par.k/2;
	//col = par.k - par.k/2;
	//c = 0;
	//for(i = 0; i < par.n; ++i){
	//	for(j = c; (j < col) && (j < par.n); ++j){
	//		printf("%lf ", A[i*par.k + (par.k/2 - i) + j]);
	//	}
	//	col++;
	//	numZeros--;
	//	if(numZeros < 0){
	//		c++;
	//	}
	//	printf("\n");
	//}
//
//	//numZeros = par.k/2;
//	//for(i = 0; i < (par.n + numZeros*2); ++i){
//	//	printf("%.1lf\n", B[i]);
	//}

	gradienteConjugado(A, B, par);

	//numZeros = par.k/2;
	//col = par.k - par.k/2;
	//c = 0;
	//for(i = 0; i < par.n; ++i){
	//	for(j = c; (j < col) && (j < par.n); ++j){
	//		printf("%lf ", A[i*par.k + (par.k/2 - i) + j]);
	//	}
	//	col++;
	//	numZeros--;
	//	if(numZeros < 0){
	//		c++;
	//	}
	//	printf("\n");
	//}
//
//	//numZeros = par.k/2;
//	//for(i = 0; i < (par.n + numZeros*2); ++i){
//	//	printf("%.1lf\n", B[i]);
	//}
	free(A);
	free(B);

	return(0);
}

