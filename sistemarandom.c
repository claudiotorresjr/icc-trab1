
/**
 * @file sistemarandom.c
 * @author GRR20176143 Cláudio Torres Júnior
 * @author GRR20171607 Gabriela Stein
 * @date 16 Sep 2018
 * @brief Arquivo de implementação das funções de geração do sistema linear.
 */

#include <stdlib.h>
#include "sistemarandom.h"

/**
 * @brief Função que gera os coeficientes de um sistema linear k-diagonal
 * @param i coordenadas do elemento a ser calculado (0<=i,j<n)
 * @param j coordenadas do elemento a ser calculado (0<=i,j<n)
 * @param k numero de diagonais da matriz A
 */

inline double generateRandomA( unsigned int i, unsigned int j, unsigned int k )
{
  static double invRandMax = 1.0 / (double)RAND_MAX;
  return ((i==j)?((double)(k<<1)):(1.0))  * ((double)rand() * invRandMax);
}

/**
 * @brief Função que gera os termos independentes de um sistema linear k-diagonal
 * @param k numero de diagonais da matriz A
 */

inline double generateRandomB( unsigned int k )
{
  static double invRandMax = 1.0 / (double)RAND_MAX;
  return ((double)(k<<2)) * ((double)rand() * invRandMax);
}