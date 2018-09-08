#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
//#include "utils.h"
#include <unistd.h> // *POSIX* Para o getopt() original
#include <getopt.h> // *GNU* Para o getopt_long()

#include "sistemarandom.h"

typedef struct parametro{
	long long int n, //dimensao
				  k, //n de diagonais
		          i; //max de iteracao
	double p, //pré-condicionador
		   e; //erro aproximado
	char *o;
}parametro;

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
	exit(0);
}

int opcoes(int argc, char *argv[], parametro *par){
	int opt;

	while((opt = getopt(argc, argv, "n:k:p:i:e:o:")) != -1){
		switch(opt){
			case 'n':
				if(isdigit(*argv[optind - 1]))
					par->n = atoi(optarg);
				else{
					printf("-n precisa de um parametro numerico.\n");
					exit(0);
				}
				if(par->n <= 10){
					printf("n precisa ser > 10.\n");
					exit(0);
				}
				break;
			case 'k':
				if(isdigit(*argv[optind - 1]))
					par->k = atoi(optarg);
				else{
					printf("-k precisa de um valor associado numerico.\n");
					exit(0);
				}
				if( (par->k > par->n) || (par->k <= 1) || ((par->k % 2) == 0) ){
					printf("k deve ser menor que n, ser > 1 e impar.\n");
					exit(0);
				}
				break;
			case 'p':
				if(isdigit(*argv[optind - 1]))
					par->p = atof(optarg);
				else{
					printf("-p precisa de um valor associado numerico.\n");
					exit(0);
				}
				if(par->p >= 0.0 && par->p > 2.0){
					printf("p precisa estar: 0 <= p <= 2.\n");
					exit(0);
				}
				break;
			case 'i':
				if(isdigit(*argv[optind - 1]))
					par->i = atoi(optarg);
				else{
					printf("-i precisa de um valor associado numerico.\n");
					exit(0);
				}
				break;
			case 'e':
				if(strcmp(optarg, "-o") == 0){
					par->e = 0.0;
					par->o = argv[optind];
				}
				else
					par->e = atof(optarg);
				break;
			case 'o':
				par->o = optarg;
				break;
			case '?':
				printf("----Parametro passado incorreto----\n");
				ajuda();
				break;
		}
	}  
	return 1; 
}

int main (int argc, char *argv[])
{
	double *A, *B;
	int sair = 0;
	parametro par;
	srand(20182);

	if(argc != 13 && argc != 12){
		ajuda();
	}
	while(sair != 1){
		sair = opcoes(argc, argv, &par);
	}
	//printf("%lld, %lld, %lf, %lld, %lf, %s\n", par.n, par.k, par.p, par.i, par.e, par.o);

	A = (double*)malloc(par.n*par.n*sizeof(double));
	B = (double*)malloc(par.n*sizeof(double));

	for(long long int i = 0; i < par.n; i++){
		for(long long int j = 0; j < par.n; j++){
			if(abs(i - j) > par.k/2){
				A[i*par.n + j] = 0;
			}else{
				A[i*par.n + j] = generateRandomA(i, j, par.k);
			}
		}
	}

	for(long long int i = 0; i < par.n; i++){
		B[i] = generateRandomB(par.k);
	}

	for(long long int i = 0; i < par.n; i++){
		for(long long int j = 0; j < par.n; j++){
			printf("%lf ", A[i*par.n + j]);
		}
		printf("= %lf", B[i]);
		printf("\n");
	}
}

