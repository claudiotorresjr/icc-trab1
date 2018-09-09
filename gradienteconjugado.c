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
					M[i*tam + j] = 1 / A[i*tam + j];
				}else{
					M[i*tam + j] = 0.0;
				}
			}
		}
	}
	//}else if(p == 1.0){//pré-condicionador de Gauss-Seidel

	//}else{//pré-condicionador SSOR

	//}
}

int gradienteConjugado(double *A, double *B, parametro par){
	double *M = (double*)malloc(par.n*par.n*sizeof(double)); //pre-condicionador
	double *X = (double*)malloc(par.n*sizeof(double)); 		//vetor de 'chutes' iniciais
	double *residuo = (double*)malloc(par.n*sizeof(double));//residuo
	double *r = (double*)malloc(par.n*sizeof(double));		
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
	multMatVet(M, B, par.n, v);

	//y = M-¹ + r
	multMatVet(M, r, par.n, y);

	//aux = y^t * r
	aux = multVetVet(y, r, par.n);

	norma[0] = sqrtf(aux); //Norma de b (erro inicial pois X = 0)
	//it == 1 pois a it 0 foi feito fora do for
	for(int it = 1; it < par.i; it++){
		//z = A*v
		multMatVet(A, v, par.n, z);

		//s = aux /v^t * z
		s = aux/multVetVet(v, z, par.n);

		//x = x + s*v
		for(i = 0; i < par.n; i++)
			X[i] = X[i] + s*v[i];

		//r = r - s*z
		for(i = 0; i < par.n; i++)
			r[i] = r[i] - s*z[i];

		//y = M-¹ * r
		multMatVet(M, r, par.n, y);

		//aux1 = r^t * r
		aux1 = multVetVet(r, r, par.n);

		//calculo do residuo -> residuo = B - A*X
		multMatVet(A, X, par.n, residuo); //residuo = A*X
		for(i = 0; i < par.n; i++){ //residuo = B - residuo 
			residuo[i] = B[i] - residuo[i];
		}

		//for(i = 0; i < par.n; i++)
		//	printf("%lf\n", residuo[i]);

		norma[it] = sqrtf(multVetVet(residuo, residuo, par.n)); //norma euclidiana do residuo
		erroAproximado = fabs((norma[it - 1] - norma[it]) / norma[it]); //erro aproximado absoluto
		if(erroAproximado < par.e){
			//achou resultado
			//encerra();
			for(i = 0; i < par.n; i++)
				printf("%lf ", X[i]);
			printf("\n ");
			printf("%d ", it);
			printf("\n ");
			return 0;
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
	return -1;
}