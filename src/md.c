#include "stdlib.h"
#include <time.h>
#include <stdio.h>
#include "dynamics.h"

#define PASO 0.02
#define PASO2 0.0002          // esto es paso al cuadrado sobre dos

// pasos a seguir:
// calcular la energia (chequear que sea cte)
// calcular la distribucion radial
// calcular presion

int main(int argc, char **argv) {

  hola();
  // Inicializamos numero de particulas
  int n_part = 27;
  float l = 3;
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

  float r_cut = 2.5;   // r_cut = (2.5 * sigma), tomamos sigma = 1

  for (int n = 0; n < n_part; n++) {
    F_tot(n, n_part, l, pos_x_ant, pos_y_ant, pos_z_ant, r_cut, fuerza_x_ant, fuerza_y_ant, fuerza_z_ant);
  }



  for (int i = 0; i < n_part; i++) {
    printf("particula = %i,   ", i);
    printf("\n");
    printf("pos x es %f,    ", pos_x_ant[i]);
    printf("vel x es %f,    ", vel_x_ant[i]);
    printf("fuerza x es %f,    ", fuerza_x_ant[i]);
    printf("\n");
    printf("pos y es %f,    ", pos_y_ant[i]);
    printf("vel y es %f,    ", vel_y_ant[i]);
    printf("fuerza y es %f,    ", fuerza_y_ant[i]);
    printf("\n");
    printf("pos z es %f,    ", pos_z_ant[i]);
    printf("vel z es %f,    ", vel_z_ant[i]);
    printf("fuerza z es %f,    ", fuerza_z_ant[i]);
    printf("\n");
    printf("\n");
  }
  return 0;
}
