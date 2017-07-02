#include "stdlib.h"
#include <time.h>
#include <stdio.h>
#include "dynamics.h"

#define PASO 0.02
#define PASO2 0.0002          // esto es paso al cuadrado sobre dos
#define R_CUT2 6.25

// pasos a seguir:
// calcular la energia (chequear que sea cte)
// calcular la distribucion radial
// calcular presion

int main(int argc, char **argv) {

  hola();
  // Inicializamos numero de particulas
  int n_part = 27;
  float l = 3;
  int steps = 1;
  float T = 0.728;
  srand(time(NULL));

  // Inicializamos arrays de energias cinetica, potencial y total
  float* potential = calloc(steps, sizeof(float));
  float* kinetic = calloc(steps, sizeof(float));
  float* E_total = calloc(steps, sizeof(float));
  float* temperature = calloc(steps, sizeof(float));

  // Inicializamos array de posiciones y velocidades
  // unos para un tiempo anterior y otros para el posterior
  float* pos_x_ant = calloc(n_part, sizeof(float));
  float* pos_y_ant = calloc(n_part, sizeof(float));
  float* pos_z_ant = calloc(n_part, sizeof(float));

  float* vel_x_ant = calloc(n_part, sizeof(float));
  float* vel_y_ant = calloc(n_part, sizeof(float));
  float* vel_z_ant = calloc(n_part, sizeof(float));

  float* fuerza_x_ant = calloc(n_part, sizeof(float));
  float* fuerza_y_ant = calloc(n_part, sizeof(float));
  float* fuerza_z_ant = calloc(n_part, sizeof(float));


  float* pos_x_post = calloc(n_part, sizeof(float));
  float* pos_y_post = calloc(n_part, sizeof(float));
  float* pos_z_post = calloc(n_part, sizeof(float));

  float* vel_x_post = calloc(n_part, sizeof(float));
  float* vel_y_post = calloc(n_part, sizeof(float));
  float* vel_z_post = calloc(n_part, sizeof(float));

  float* fuerza_x_post = calloc(n_part, sizeof(float));
  float* fuerza_y_post = calloc(n_part, sizeof(float));
  float* fuerza_z_post = calloc(n_part, sizeof(float));


  initizalize_pos(n_part, l, pos_x_ant, pos_y_ant, pos_z_ant);
  initizalize_vel(n_part, vel_x_ant, vel_y_ant, vel_z_ant);


  for (int n = 0; n < n_part; n++) {
    F_tot(n, n_part, l, pos_x_ant, pos_y_ant, pos_z_ant, R_CUT2, fuerza_x_ant, fuerza_y_ant, fuerza_z_ant);
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

  kinetic_temperature(n_part, 0, vel_x_ant, vel_y_ant, vel_z_ant, temperature, kinetic);
  potential_energy(n_part, 0, l, pos_x_ant, pos_y_ant, pos_z_ant, 6.25, potential);
  total_energy(0, kinetic, potential, E_total);

  printf("ENERGIAS\n");
  printf("cinetica = %f, potencial = %f, total = %f\n", kinetic[0], potential[0], E_total[0]);
  printf("temperatura = %f\n", temperature[0]);
  printf("\n");
  printf("\n");
  return 0;
}
