#include "stdlib.h"
#include <time.h>
#include <stdio.h>
#include "dynamics.h"

int main(int argc, char **argv) {

  hola();
  // Inicializamos numero de particulas
  int n_part = 8;
  float l = 2;
  // Inicializamos array de posiciones y velocidades
  // unos para un tiempo anterior y otros para el posterior
  float* pos_x_ant = malloc(n_part * sizeof(float));
  float* pos_y_ant = malloc(n_part * sizeof(float));
  float* pos_z_ant = malloc(n_part * sizeof(float));

  float* vel_x_ant = malloc(n_part * sizeof(float));
  float* vel_y_ant = malloc(n_part * sizeof(float));
  float* vel_z_ant = malloc(n_part * sizeof(float));

  float* fuerza_x_ant = malloc(n_part * sizeof(float));
  float* fuerza_y_ant = malloc(n_part * sizeof(float));
  float* fuerza_z_ant = malloc(n_part * sizeof(float));


  float* pos_x_post = malloc(n_part * sizeof(float));
  float* pos_y_post = malloc(n_part * sizeof(float));
  float* pos_z_post = malloc(n_part * sizeof(float));

  float* vel_x_post = malloc(n_part * sizeof(float));
  float* vel_y_post = malloc(n_part * sizeof(float));
  float* vel_z_post = malloc(n_part * sizeof(float));

  float* fuerza_x_post = malloc(n_part * sizeof(float));
  float* fuerza_y_post = malloc(n_part * sizeof(float));
  float* fuerza_z_post = malloc(n_part * sizeof(float));


  initizalize_pos(n_part, l, pos_x_ant, pos_y_ant, pos_z_ant);
  initizalize_vel(n_part, vel_x_ant, vel_y_ant, vel_z_ant);

  for (int i = 0; i < n_part; i++) {
    printf("pos x es %f\n", pos_x_ant[i]);
    printf("vel x es %f\n", vel_x_ant[i]);
  }

  return 0;
}
