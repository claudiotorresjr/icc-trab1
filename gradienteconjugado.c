#include "gradienteconjugado.h"

void multMatMat(double *pri, double *sec, long int tam, double *mult){
	long int i, j, k;
	double soma = 0.0;

	for (i = 0; i < tam; i++){
		for (j = 0; j < tam; j++){
			for (k = 0; k < tam; k++){
				soma = soma + pri[i*tam + k]*sec[k*tam + j];
			}
			mult[i*tam + j] = soma;
			soma = 0.0;
		}
	}
}

void multMatVet(double *pri, double *sec, long int tam, double *mult){
	long int i, j, k;
	double soma = 0.0;

	for (i = 0; i < tam; i++){
		for (j = 0; j < tam; j++){
			for (k = 0; k < tam; k++){
				soma = soma + pri[i*tam + k]*sec[k];
			}
			mult[i] = soma;
			soma = 0.0;
		}
	}
}

double multVetVet(double *pri, double *sec, long int tam){
	long int i, j;
	double soma = 0.0;

	for (i = 0; i < tam; i++){
		soma = soma + pri[i]*sec[i];
	}
	return soma;
}

void trasformaSistema(double *A, double *B, long int tam){
	double *T = (double*)malloc(tam*tam*sizeof(double));
	double *mult = (double*)malloc(tam*tam*sizeof(double));
	double *multv = (double*)malloc(tam*sizeof(double));

	long int i, j;

	transposta(A, T, tam);
	multMatMat(T, A, tam, mult);
	multMatVet(T, B, tam, multv);
	for(i = 0; i < tam; i++){
		for(j = 0; j < tam; j++){
			A[i*tam + j] = mult[i*tam + j];
			B[i] = multv[i];
		}
	}
}

void transposta(double *A, double *T, long int tam){
	long int i, j;

	for(i = 0; i < tam; i++){
		for(j = 0; j < tam; j++){
			T[j*tam + i] = A[i*tam + j];
		}
	}
}

void preCondicionador(double p, double *M, double *A, long int tam){
	long int i, j;
	double *L = (double*)malloc(tam*tam*sizeof(double)); //matriz lower
	double *U = (double*)malloc(tam*tam*sizeof(double)); //matriz upper
	double *D = (double*)malloc(tam*tam*sizeof(double)); //matriz diagonal
	double *R1 = (double*)malloc(tam*tam*sizeof(double)); //matriz resposta1
	double *R2 = (double*)malloc(tam*tam*sizeof(double)); //matriz resposta2

	if(p == 0.0){//sem pre-condicionador
		for(i = 0; i < tam; i++){
			for(j = 0; j < tam; j++){
				if(i == j){
					M[i*tam + j] = 1.0;
				}else{
					M[i*tam + j] = 0.0;
				}
			}
		}
	}else if(p > 0.0 && p < 1.0){//pre-condicionador de jacobi
		for(i = 0; i < tam; i++){
			for(j = 0; j < tam; j++){
				if(i == j){
					M[i*tam + j] = A[i*tam + j];
				}else{
					M[i*tam + j] = 0.0;
				}
			}
		}
	}else if(p == 1.0){//pré-condicionador de Gauss-Seidel

		criaMatrizes(A, L, D, U, tam);
		for(i = 0; i < tam; i++){
			for(j = 0; j < tam; j++){
					R1[i*tam + j] = D[i*tam + j] + L[i*tam + j];
			}
		}
		for(i = 0; i < tam; i++){
			for(j = 0; j < tam; j++){
					R2[i*tam + j] = D[i*tam + j] + U[i*tam + j];
			}
		}
		for(i = 0; i < tam; i++){
			for(j = 0; j < tam; j++){
				if(i == j){
					D[i*tam + j] = 1/A[i*tam + j];
				}else{
					D[i*tam + j] = 0.0;
				}
			}
		}
		multMatMat(R1, D, tam, U); //matriz U sendo usada como intermediaria da operação
		multMatMat(U, R2, tam, M);

	}else{//pré-condicionador SSOR
		criaMatrizes(A, L, D, U, tam);
		for(i = 0; i < tam; i++){
			for(j = 0; j < tam; j++){
					R1[i*tam + j] = D[i*tam + j] + p*L[i*tam + j];
			}
		}
		for(i = 0; i < tam; i++){
			for(j = 0; j < tam; j++){
					R2[i*tam + j] = D[i*tam + j] + p*U[i*tam + j];
			}
		}
		for(i = 0; i < tam; i++){
			for(j = 0; j < tam; j++){
				if(i == j){
					D[i*tam + j] = 1/A[i*tam + j];
				}else{
					D[i*tam + j] = 0.0;
				}
			}
		}
		multMatMat(R1, D, tam, U); //matriz U sendo usada como intermediaria da operação
		multMatMat(U, R2, tam, M);
	}
}

double maxVetor(double *V, long int tam)
{
	double maximo = V[0];
	for (int i = 0; i < tam; i++)
	{
		if (V[i] > maximo)
			maximo = V[i];
	}
	return maximo;
}

/*void imprime_dados()
{
	# login1 Nome1
	# login2 Nome2
	#
	# iter 1: <||x||>
	# iter 2: <||x||>
	# ...
	# iter k: <||x||>
	# residuo: <||r||>
	# Tempo PC: <tempo para cálculo do pré-condicionador>
	# Tempo iter: <tempo para resolver uma iteração do método>
	# Tempo residuo: <tempo para calcular o residuo do SL>
	#
	n
	x_1 x_12 ... x_n
} */

void Cholesky(double *M, double *a, double *b)
{
	
}

void criaMatrizes(double *A, double *L, double *U, double *D, long int tam)
{
	long int i, j;
	//cria a matriz upper
		for(i = 0; i < tam; i++){
			for(j = 0; j < tam; j++){
				if(j > i){
					U[i*tam + j] = A[i*tam + j];
				}else{
					U[i*tam + j] = 0.0;
				}
			}
		}
		//cria a matriz lower
		for(i = 0; i < tam; i++){
			for(j = 0; j < tam; j++){
				if(j < i){
					L[i*tam + j] = A[i*tam + j];
				}else{
					L[i*tam + j] = 0.0;
				}
			}
		}
		//cria a matriz diagonal
		for(i = 0; i < tam; i++){
			for(j = 0; j < tam; j++){
				if(j == i){
					D[i*tam + j] = A[i*tam + j];
				}else{
					D[i*tam + j] = 0.0;
				}
			}
		}
}

int gradienteConjugado(double *A, double *B, parametro par){
	int convergiu = 0;
	double *M = (double*)malloc(par.n*par.n*sizeof(double)); //pre-condicionador
	double *X = (double*)malloc(par.n*sizeof(double)); 		//vetor de 'chutes' iniciais
	double *Xant = (double*)malloc(par.n*sizeof(double));
	double *r = (double*)malloc(par.n*sizeof(double));		//residuo
	double *v = (double*)malloc(par.n*sizeof(double));		
	double *z = (double*)malloc(par.n*sizeof(double));		
	double *y = (double*)malloc(par.n*sizeof(double));	
	double *T = (double*)malloc(par.n*sizeof(double));		//usado para calcular as transpostas
	double *norma = malloc(par.n*sizeof(double));           //vetor com os erros (usado para tirar o max erro)

	double s, aux, aux1, m, erroAproximado;

	long int i, j;

	preCondicionador(par.p, M, A, par.n);

	//X0 = 0
	for(i = 0; i < par.n; i++)
		X[i] = 0;

	//r = B
	for(i = 0; i < par.n; i++)
		r[i] = B[i];

	//v = M-¹ * B
	for(i = 0; i < par.n; i++)
		v[i] = B[i] / M[i*par.n + i];


	//y = M-¹ * r
	for(i = 0; i < par.n; i++)
		y[i] = r[i] / M[i*par.n + i];

	//aux = y^t * r
	aux = multVetVet(y, r, par.n);

	norma[0] = sqrtf(aux); //Norma de b (erro inicial pois X = 0)

	//it == 1 pois a it 0 foi feito fora do for
	for(int it = 1; it < par.i; it++){
		//z = A*v
		multMatVet(A, v, par.n, z);

		//s = aux /v^t * z
		s = aux/multVetVet(v, z, par.n);

		//salva o vetor
		for(i = 0; i < par.n; i++)
			Xant[i] = X[i];

		//x = x + s*v
		for(i = 0; i < par.n; i++)
			X[i] = X[i] + s*v[i];

		//r = r - s*z
		for(i = 0; i < par.n; i++)
			r[i] = r[i] - s*z[i];  //calculo do residuo

		//y = M-¹ * r
		for(i = 0; i < par.n; i++)
		y[i] = r[i] / M[i*par.n + i];

		//aux1 = r^t * r
		aux1 = multVetVet(r, r, par.n);

		//calculo do maior elemento do vetor X e do vetor X anterior

		erroAproximado = fabs((maxVetor(X, par.n) - maxVetor(Xant, par.n)) / maxVetor(X, par.n)); //erro aproximado absoluto

		//ESSA PARTE AQUI

		//calculo do residuo -> residuo = B - A*X
		/*multMatVet(A, X, par.n, residuo); //residuo = A*X
		for(i = 0; i < par.n; i++){ //residuo = B - residuo 
			residuo[i] = B[i] - residuo[i];
		}*/

		//for(i = 0; i < par.n; i++)
		//	printf("%lf\n", residuo[i]);

		norma[it] = sqrtf(multVetVet(r, r, par.n)); //norma euclidiana do residuo


		//ATE AQUI

		//printf("%lf\n", erroAproximado);

		if(erroAproximado < par.e && !convergiu){
			//achou resultado
			//encerra();
			for(i = 0; i < par.n; i++)
				printf("%lf ", X[i]);
			printf("\n ");
			printf("%d ", it);
			printf("\n ");
			if (par.op == 0) //se falso, para as iterações
				return 0;
			else
				convergiu = 1;
		}
		//aux1 = y^t * r;
		aux1 = multVetVet(y, r, par.n);

		//m = aux1/aux
		m = aux1/aux;

		//aux = aux1
		aux = aux1;

		//v = y + m*v
		for(i = 0; i < par.n; i++)
			v[i] = y[i] + m*v[i];
	}
	if (!convergiu)
		fprintf(stderr, "O método não convergiu!\n");
	return -1;
}