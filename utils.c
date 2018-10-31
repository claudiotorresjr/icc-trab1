
/**
 * @file utils.c
 * @author GRR20176143 Cláudio Torres Júnior
 * @author GRR20171607 Gabriela Stein
 * @date 16 Sep 2018
 * @brief Arquivo de implementação da função para medição do tempo.
 */

#include "utils.h"

// Retorna tempo em milisegundos
double timestamp(void)
{
  struct timeval tp;
  gettimeofday(&tp, NULL);
  return((double)((tp.tv_sec*1000.0 + tp.tv_usec/1000.0)*0.001));
}

