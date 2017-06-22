#include "stdlib.h"
#include <time.h>
#include <stdio.h>
#include "dynamics.h"


int main(int argc, char **argv) {

  hola();
  // Inicializamos numero de particulas
  int n_part = 1;

  // Inicializamos array de posiciones y velocidades
  // unos para un tiempo anterior y otros para el posterior
  float* pos_x_ant = malloc(n_part * sizeof(float));
  float* pos_y_ant = malloc(n_part * sizeof(float));
  float* pos_z_ant = malloc(n_part * sizeof(float));

  float* vel_x_ant = malloc(n_part * sizeof(float));
  float* vel_y_ant = malloc(n_part * sizeof(float));
  float* vel_z_ant = malloc(n_part * sizeof(float));

  float* pos_x_post = malloc(n_part * sizeof(float));
  float* pos_y_post = malloc(n_part * sizeof(float));
  float* pos_z_post = malloc(n_part * sizeof(float));

  float* vel_x_post = malloc(n_part * sizeof(float));
  float* vel_y_post = malloc(n_part * sizeof(float));
  float* vel_z_post = malloc(n_part * sizeof(float));



  return 0;
}
