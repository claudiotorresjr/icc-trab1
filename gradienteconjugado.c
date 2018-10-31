
/**
 * @file gradienteconjugado.c
 * @author GRR20176143 Cláudio Torres Júnior
 * @author GRR20171607 Gabriela Stein
 * @date 16 Sep 2018
 * @brief Arquivo de implementação das principais funções.
 */

#include "gradienteconjugado.h"
#include "utils.h"

/**
 * @brief Função que multiplica duas matrizes
 * @param pri Primeiro parametro: Matriz
 * @param sec Segundo parametro: Matriz
 * @param mult Resultado da multiplicação
 * @param tam Ordem da Matriz
*/

void multMatMat(double *pri, double *sec, long int dgn, long int tam, double *mult){
	long int i, j, k;
	double soma = 0.0;

	long int numZerosJ = dgn/2;	
	long int col = dgn - dgn/2;
	long int cJ = 0;

	long int cK = 0;
	long int numZerosK = dgn/2;	
	long int colK;
	//long int c;

	for(i = 0; i < tam; ++i){
		numZerosK = dgn/2;
		colK = dgn - dgn/2;
		//if(cK > cJ){
		//	cK = cJ;
		//}
		for(j = 0; j < tam; ++j){
			for(k = cK; (k < col) && (k < colK); ++k){
				soma = soma + pri[i*dgn + (dgn/2 - i) + k]*sec[k*dgn + (dgn/2 - k) + j];
			}
			mult[i*(dgn*2 - 1) + ((dgn*2 - 1)/2 - i) + j] = soma;
			soma = 0.0;
			colK++;
			numZerosK--;
			if(numZerosK < 0){
				if(abs(numZerosK) > abs(numZerosJ)){
					cK = abs(numZerosK);
				}
			}
		}
		col++;
		numZerosJ--;
		if(numZerosJ < 0){
			cJ++;
			cK = cJ;
		}else{
			cK = 0;
		}
	}
}

/**
 * @brief Função que multiplica uma matriz por um vetor
 * @param pri Primeiro parametro: Matriz
 * @param sec Segundo parametro: Vetor
 * @param mult Resultado da multiplicação
 * @param tam Ordem da Matriz
*/

void multMatVet(double *pri, double *sec, long int inicio, long int dgn, long int tam, double *mult){
	long int i, j, c = 0;
	long int numZeros = dgn/2;
	long int col = dgn - dgn/2;	
	double soma = 0.0;

	for(i = 0; i < tam; i++){
		for(j = c; (j < col) && (j < tam); j++){
			soma = soma + pri[i*dgn + (dgn/2 - i) + j]*sec[inicio + j];
		}
		mult[inicio + i] = soma;
		soma = 0.0;
		col++;
		numZeros--;
		if(numZeros < 0){
			c++;
		}
	}
}

/**
 * @brief Função que multiplica dois vetores
 * @param pri Primeiro parametro: Vetor
 * @param sec Segundo parametro: Vetor
 * @param tam Ordem da Matriz
*/

double multVetVet(double *pri, double *sec, long int dgn, long int tam){
	long int i;
	double soma = 0.0;
	long int numZeros = dgn/2;

	for (i = 0; i < (tam + numZeros); i++){
		soma = soma + pri[i]*sec[i];
	}
	return soma;
}

/**
 * @brief Função que transforma um sistema em um sistema simetrico e positivo definido
 * @param A matriz original
 * @param B vetor de termos independentes
 * @param tam Ordem da Matriz
*/

void trasformaSistema(double *A, double *B, double *Atf, double *Btf, parametro par){
	long int tam = par.n*par.k;

	double *T = (double*)malloc(tam*sizeof(double));

	memset(T, 0, sizeof(*T));

	transposta(A, T, par);
	multMatMat(T, A, par.k, par.n, Atf);
	multMatVet(T, B, par.k/2, par.k, par.n, Btf);

	free(T);
}

/**
 * @brief Função que calcula a transposta de uma matriz
 * @param A matriz original
 * @param T matriz transposta
 * @param tam Ordem da Matriz
*/

void transposta(double *A, double *T, parametro par){
	long int i, j;
	long int numZeros = par.k/2;
	long int col = par.k - par.k/2;
	long int c = 0;

	for(i = 0; i < par.n; i++){
		for(j = c; (j < col) && (j < par.n); j++){
			T[j*par.k + (par.k/2 - j) + i] = A[i*par.k + (par.k/2 - i) + j];
		}
		col++;
		numZeros--;
		if(numZeros < 0){
			c++;
		}
	}
}

/**
 * @brief Função que condiciona uma matriz
 * @param A matriz original
 * @param M matriz resultante do uso do pré condicionador
 * @param p Indica o pré-condicionador a ser utilizado
 * @param tam Ordem da Matriz
*/

/*void preCondicionador(double p, double *M, double *A, long int tam){
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
}*/

/**
 * @brief Função que encontra o maior valor de um vetor
 * @param V Vetor qualquer
 * @param tam Ordem da Matriz
 * @return Retorna o maior valor
*/

double maxVetor(double *V, parametro par)
{	
	long int i;
	long int numZeros = par.k/2;
	double maximo = V[0];
	for (i = 0; i < (par.n + numZeros); ++i)
	{
		if (V[i] > maximo)
			maximo = V[i];
	}
	return maximo;
}

/**
 * @brief Função que imprime a saída de dados
 * @param erroIT Vetor que contem o erro em todas iterações
 * @param X Vetor solução
 * @param norma Norma do resíduo
 * @param pc Tempo para cálculo do pré-condicionador
 * @param it Tempo para resolver uma iteração
 * @param r Tempo usado para calcular o residuo do SL
 * @param par Struct contendo os parametros passados na execução do programa
 * @param iter Quantidade de iterações usadas pelo método
*/

void imprime_dados(double *erroIt, double *X, double norma, double pc, double it, double r, parametro par, long int iter)
{	
	long int numZeros = par.k/2;
	FILE *arqOut;
	arqOut = fopen(par.o, "w");
	if(arqOut == NULL){
		fprintf(stderr, "erro ao criar arquivo\n");
		exit(1);
	}

	fprintf(arqOut, "# ctj17 Cláudio Torres Júnior\n"); //# login1 Nome1
	fprintf(arqOut, "# gs17 Gabriela Stein\n#\n"); //# login2 Nome2
	for(long int i = 1; i <= iter; i++){
		fprintf(arqOut, "# iter %ld: <||%lf||>\n", i, erroIt[i]); //# iter k: <||x||>
	}
	fprintf(arqOut, "# residuo: <||%lf||>\n", norma); //# residuo: <||r||>
	fprintf(arqOut, "# Tempo PC: <%lf>\n", pc); //# Tempo PC: <tempo para cálculo do pré-condicionador>
	fprintf(arqOut, "# Tempo iter: <%lf>\n", it); //# Tempo iter: <tempo para resolver uma iteração do método>
	fprintf(arqOut, "# Tempo residuo: <%lf>\n#\n%ld\n", r, par.n); //# Tempo residuo: <tempo para calcular o residuo do SL> 
	for(long int i = par.k/2; i < (par.n + numZeros); i++){
		fprintf(arqOut, "%.15g ", X[i]); //x_1 x_12 ... x_n
	}
	fprintf(arqOut, "\n");
} 

/**
 * @brief Função que calcula a Fatoração Incompleta de Cholesky
 * @param M Matriz condionada usando Gauss Seidel ou SSOR
 * @param tam Ordem da Matriz
*/

/*void Cholesky(double *M, long int tam)
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
}*/

/**
 * @brief Função que aloca e atribui as matrizes criadas para os pré condionantes Gauss Seidel e SSOR
 * @param A Matriz original
 * @param L Matriz inferior (lower)
 * @param U Matriz superior (upper)
 * @param D Matriz diagonal
 * @param tam Ordem da Matriz
*/

/*void criaMatrizes(double *A, double *L, double *U, double *D, long int tam)
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

				if(j < i){
					L[i*tam + j] = A[i*tam + j];
				}else{
					L[i*tam + j] = 0.0;
				}

				if(j == i){
					D[i*tam + j] = A[i*tam + j];
				}else{
					D[i*tam + j] = 0.0;
				}
			}
		}
}*/

/**
 * @brief Função que libera os vetores criados durante a execução do programa
 * @param M Matriz condicionada
 * @param X Vetor solução
 * @param Xant Vetor solução da iteração anterior
 * @param r Resíduo
 * @param v M^(-1)*B
 * @param z Resultado para a matriz verdadeira * vetor resposta de M
 * @param y vetor soluçao para o sistema M*y = r
 * @param T matriz transposta
 * @param erroAproximadoR vetor para o erro aproximado relativo
 * @param erroAproximadoA vetor para o erro aproximado absoluto
 * @param erroIt vetor para o erro maximo de cada iteração
*/

void liberaVet(double *M, double *X, double *r, double *v, double *z, 
	double *y, double *Xant, double *erroAproximadoA, double *erroIt){

	free(M);
	printf("M\n");
	free(X);
	printf("X\n");
	free(r);
	printf("r\n");
	free(v);
	printf("v\n");
	free(z);
	printf("z\n");
	free(y);
	printf("y\n");
	free(Xant);
	printf("Xant\n");
	free(erroAproximadoA);
	printf("erroAproximadoA\n");
	free(erroIt);
	printf("erroIt\n");
}

/**
 * @brief Função que calcula a solução do sistema linear pelo método do Gradiente Conjugado
 * @param A Matriz do sistema 
 * @param B Vetor solução
 * @param par Struct contendo os parametros passados durante a chamada do programa
*/

int gradienteConjugado(double *A, double *B, parametro par){
	int convergiu = 0;
	long int numZeros = par.k/2;

	double *M = (double*)malloc((par.n + numZeros)*sizeof(double)); //pre-condicionador
	double *X = (double*)malloc((par.n + numZeros)*sizeof(double)); 		//vetor de 'chutes' iniciais
	double *Xant = (double*)malloc((par.n + numZeros)*sizeof(double));
	double *r = (double*)malloc((par.n + numZeros)*sizeof(double));		//residuo
	double *v = (double*)malloc((par.n + numZeros)*sizeof(double));		
	double *z = (double*)malloc((par.n + numZeros)*sizeof(double));		
	double *y = (double*)malloc((par.n + numZeros)*sizeof(double));	
	////double *erroAproximadoR = malloc(par.n*sizeof(double));           //vetor com o erro relativo
	double *erroAproximadoA = malloc((par.n + numZeros)*sizeof(double));	//vetor com o erro absoluto
	double *erroIt = malloc(par.i*sizeof(double));	//vetor com o erro max absoluto em cada iteraçao
	double *Atf = (double*)malloc((par.k*2 - 1)*par.n*sizeof(double)); 
	double *Btf = (double*)malloc((par.n + numZeros)*sizeof(double)); 		

	double s, aux, aux1, m, norma; //Xprox;

	long int i, j, it;

	memset(Atf, 0, sizeof(*Atf));
	memset(M, 0, sizeof(*M));

	tempo t_pc; //struct para usar o timestamp precodicionador
	tempo t_it; //struct para usar o timestamp cada iteraçao
	tempo t_r; //struct para usar o timestamp residuo

	t_pc.ini = timestamp();
	//transforma A em uma matriz positiva simetrica
	trasformaSistema(A, B, Atf, Btf, par);
	//acha pre-condicionador
	//preCondicionador(par.p, M, Atf, par.n);
	for(i = 0; i < (par.n + numZeros); ++i){
		M[i + par.k/2] = Atf[i*(par.k*2 - 1) + ((par.k*2 - 1)/2 - i) + i];
	}

	t_pc.fim = timestamp();
	t_pc.dif = t_pc.fim - t_pc.ini;

	
	for(i = 0; i < (par.n + numZeros); ++i){
		X[i] = 0; //X0 = 0
		r[i] = B[i]; //r = B
	}
		

	if(par.p < 1){
//		//v = M-¹ * B e y = M-¹ * r pois r == B
		for(i = par.k/2; i < (par.n + numZeros); ++i){
			v[i] = Btf[i] / M[i];
			y[i] = v[i];
		}
	}else{
		double soma = 0.0;
		v[0] = Btf[0] / M[0*par.n + 0];
		y[0] = v[0];
		for(i = 1; i < par.n; i++){
			soma = Btf[i];
			for(j = i - 1; j >= 0; j--){
				soma = soma - M[i*par.n + j]*v[j];
			}
			v[i] = soma / M[i*par.n + i];
			y[i] = v[i];
		}
	}

//	//aux = y^t * r
	aux = multVetVet(y, r, par.k, par.n);

//	//it == 1 pois a it 0 foi feito fora do for
	for(it = 1; it < par.i; it++){
		t_it.ini = timestamp();

//		//z = A*v
		multMatVet(Atf, v, par.k/2, (par.k*2 - 1), par.n, z);

		//numZeros = (par.k*2 - 1)/2;
		//long int col = (par.k*2 - 1) - (par.k*2 - 1)/2;
		//long int c = 0;
		//for(i = 0; i < par.n; ++i){
		//	for(j = c; (j < col) && (j < par.n); ++j){
		//		printf("%lf ", Atf[i*(par.k*2 - 1) + ((par.k*2 - 1)/2 - i) + j]);
		//	}
		//	col++;
		//	numZeros--;
		//	if(numZeros < 0){
		//		c++;
		//	}
		//	printf("\n");
		//}
//		//s = aux /v^t * z
		s = aux/multVetVet(v, z, par.k, par.n);

//		//salva o vetor
		for(i = 0; i < (par.n + numZeros); ++i)
			Xant[i] = X[i];

//		//x = x + s*v
		for(i = 0; i < (par.n + numZeros); ++i){
			X[i] = X[i] + s*v[i];
		}
//
//		//r = r - s*z
		for(i = 0; i < (par.n + numZeros); ++i){
			r[i] = r[i] - s*z[i];  //calculo do residuo
		}


//		//y = M-¹ * r
		if(par.p < 1){
			for(i = par.k/2; i < (par.n + numZeros); ++i){
				y[i] = r[i] / M[i];
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

//		//aux1 = r^t * r
//		//aux1 = multVetVet(r, r, par.n);

//		//erro aproximado absoluto
		for(i = 0; i < (par.n + numZeros); ++i){
			erroAproximadoA[i] = fabs(X[i] - Xant[i]); 
//		//	erroAproximadoR[i] = erroAproximadoA[i] / X[i];
		}
		erroIt[it] = maxVetor(erroAproximadoA, par);

		//for(i = 0; i < par.n; i++)
		//	printf("%lf\n", residuo[i]);

//		//printf("%lf\n", erroAproximado);

		if(maxVetor(erroAproximadoA, par) < par.e && !convergiu){
//			//achou resultado
			t_r.ini = timestamp();
//			//calcula a norma
			norma = sqrtf(multVetVet(r, r, par.k, par.n)); //norma euclidiana do residuo
			t_r.fim = timestamp();
			t_r.dif = t_r.fim - t_r.ini;
//			//imprime dados no arquivo
			imprime_dados(erroIt, X, norma, t_pc.dif, t_it.dif, t_r.dif, par, it);
			/*for(i = 0; i < par.n; i++)
				printf("%lf ", X[i]);
			printf("\n ");*/
			if (par.op == 0){ //se falso, para as iterações
				//liberaVet(M, X, r, v, z, y, Xant, erroAproximadoA, erroIt); 
				return 0;
			}else{
				convergiu = 1;
			}
		}

//		//aux1 = y^t * r;
		aux1 = multVetVet(y, r, par.k, par.n);

//		//m = aux1/aux
		m = aux1/aux;

//		//aux = aux1
		aux = aux1;

//		//v = y + m*v
		for(i = 0; i < (par.n + numZeros); ++i){
			v[i] = y[i] + m*v[i];
		}

	t_it.fim = timestamp();
	t_it.dif = t_it.fim - t_it.ini;
	}
	if (!convergiu){
		imprime_dados(erroIt, X, norma, t_pc.dif, t_it.dif, t_r.dif, par, it);
		fprintf(stderr, "O método não convergiu!\n");
	}

	//liberaVet(M, X, r, v, z, y, Xant, erroAproximadoA, erroIt); 
	return -1;
}