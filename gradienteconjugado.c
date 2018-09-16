#include "gradienteconjugado.h"
#include "utils.h"

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
	free(T);
	free(mult);
	free(multv);
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
	}else{//pré-condicionador de Gauss-Seidel p/ x = 1 e pré-condicionador SSOR p/ 1 < x < 2 

		criaMatrizes(A, L, U, D, tam);
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
					D[i*tam + j] = 1 / D[i*tam + j];
				}
			}
		}
		multMatMat(R1, D, tam, U); //matriz U sendo usada como intermediaria da operação
		multMatMat(U, R2, tam, M);
		Cholesky(M, tam);
	}
	free(L);
	free(U);
	free(D);
	free(R1);
	free(R2);
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

void imprime_dados(double *erroIt, double *X, double norma, double pc, double it, double r, parametro par, long int iter)
{	
	FILE *arqOut;
	arqOut = fopen(par.o, "w");
	if(arqOut == NULL){
		fprintf(stderr, "erro ao criar arquivo\n");
		exit(1);
	}

	fprintf(arqOut, "# ctj17 Cláudio Torres Júnior\n"); //# login1 Nome1
	fprintf(arqOut, "# ctj17 Cláudio Torres Júnior\n#\n"); //# login2 Nome2
	for(long int i = 1; i <= iter; i++){
		fprintf(arqOut, "# iter %ld: <||%lf||>\n", i, erroIt[i]); //# iter k: <||x||>
	}
	fprintf(arqOut, "# residuo: <||%lf||>\n", norma); //# residuo: <||r||>
	fprintf(arqOut, "# Tempo PC: <%lf>\n", pc); //# Tempo PC: <tempo para cálculo do pré-condicionador>
	fprintf(arqOut, "# Tempo iter: <%lf>\n", it); //# Tempo iter: <tempo para resolver uma iteração do método>
	fprintf(arqOut, "# Tempo residuo: <%lf>\n#\n%ld", r, par.n); //# Tempo residuo: <tempo para calcular o residuo do SL> 
	for(long int i = 0; i < par.n; i++){
		fprintf(arqOut, "%lf ", X[i]); //x_1 x_12 ... x_n
	}
} 

void Cholesky(double *M, long int tam)
{	
	long int i, j, k;
	
	for(k = 0; k < tam; k++){
		M[k*tam + k] = sqrtf(M[k*tam + k]);
		for(i = k + 1; i < tam; i++){
			if(M[i*tam + k] != 0.0){
				M[i*tam + k] = M[i*tam + k] / M[k*tam + k];
			}
		}
		for(j = k + 1; j < tam; j++){
			for(i = j; i < tam; i++){
				if(M[i*tam + j] != 0.0){
					M[i*tam + j] = M[i*tam + j] - M[i*tam + k]*M[j*tam + k];
				}
			}
		}
	}
	for(i = 0; i < tam; i++){
		for(j = i + 1; j < tam; j++){
			M[i*tam + j] = 0.0;
		}
	}	
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

void liberaVet(double *M, double *X, double *Xant, double *r, double *v, double *z, 
	double *y, double *T, double *erroAproximadoR , double *erroAproximadoA, double *erroIt){

	free(M);
	free(X);
	free(Xant);
	free(r);
	free(v);
	free(z);
	free(y);
	free(T);
	free(erroAproximadoR);
	free(erroAproximadoA);
	free(erroIt);
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
	double *erroAproximadoR = malloc(par.n*sizeof(double));           //vetor com o erro relativo
	double *erroAproximadoA = malloc(par.n*sizeof(double));	//vetor com o erro absoluto
	double *erroIt = malloc(par.i*sizeof(double));	//vetor com o erro max absoluto em cada iteraçao

	double s, aux, aux1, m, norma;

	long int i, j;

	tempo t_pc; //struct para usar o timestamp precodicionador
	tempo t_it; //struct para usar o timestamp cada iteraçao
	tempo t_r; //struct para usar o timestamp residuo

	t_pc.ini = timestamp();
	//transforma A em uma matriz positiva simetrica
	trasformaSistema(A, B, par.n);
	//acha pre-condicionador
	preCondicionador(par.p, M, A, par.n);
	t_pc.fim = timestamp();
	t_pc.dif = t_pc.fim - t_pc.ini;

	//X0 = 0
	for(i = 0; i < par.n; i++)
		X[i] = 0;

	//r = B
	for(i = 0; i < par.n; i++)
		r[i] = B[i];

	if(par.p < 1){
		//v = M-¹ * B e y = M-¹ * r pois r == B
		for(i = 0; i < par.n; i++){
			v[i] = B[i] / M[i*par.n + i];
			y[i] = v[i];
		}
	}else{
		double soma = 0.0;
		v[0] = B[0] / M[0*par.n + 0];
		y[0] = v[0];
		for(i = 1; i < par.n; i++){
			soma = B[i];
			for(j = i - 1; j >= 0; j--){
				soma = soma - M[i*par.n + j]*v[j];
			}
			v[i] = soma / M[i*par.n + i];
			y[i] = v[i];
		}
	}



	//aux = y^t * r
	aux = multVetVet(y, r, par.n);

	//it == 1 pois a it 0 foi feito fora do for
	for(int it = 1; it < par.i; it++){
		t_it.ini = timestamp();

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
		if(par.p < 1){
			for(i = 0; i < par.n; i++){
				y[i] = r[i] / M[i*par.n + i];
			}
		}else{
			double soma = 0.0;
			y[0] = r[0] / M[0*par.n + 0];
			for(i = 1; i < par.n; i++){
				soma = r[i];
				for(j = i - 1; j >= 0; j--){
					soma = soma - M[i*par.n + j]*y[j];
				}
				y[i] = soma / M[i*par.n + i];
			}
		}

		//aux1 = r^t * r
		aux1 = multVetVet(r, r, par.n);

		//erro aproximado absoluto
		for(i = 0; i < par.n; i++){
			erroAproximadoA[i] = fabs(X[i] - Xant[i]); 
			erroAproximadoR[i] = erroAproximadoA[i] / X[i];
		}
		erroIt[it] = maxVetor(erroAproximadoA, par.n);

		//for(i = 0; i < par.n; i++)
		//	printf("%lf\n", residuo[i]);

		//printf("%lf\n", erroAproximado);

		if(maxVetor(erroAproximadoR, par.n) < par.e && !convergiu){
			//achou resultado
			t_r.ini = timestamp();
			//calcula a norma
			norma = sqrtf(multVetVet(r, r, par.n)); //norma euclidiana do residuo
			t_r.fim = timestamp();
			t_r.dif = t_r.fim - t_r.ini;
			//imprime dados no arquivo
			imprime_dados(erroIt, X, norma, t_pc.dif, t_it.dif, t_r.dif, par, it);
			/*for(i = 0; i < par.n; i++)
				printf("%lf ", X[i]);
			printf("\n ");*/
			printf("%d ", it);
			printf("\n ");
			if (par.op == 0){ //se falso, para as iterações
				liberaVet(M, X, Xant, r, v, z, y, T, erroAproximadoR , erroAproximadoA, erroIt);
				return 0;
			}else 
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

	t_it.fim = timestamp();
	t_it.dif = t_it.fim - t_it.ini;
	}
	if (!convergiu)
		fprintf(stderr, "O método não convergiu!\n");

	liberaVet(M, X, Xant, r, v, z, y, T, erroAproximadoR , erroAproximadoA, erroIt);
	return -1;
}