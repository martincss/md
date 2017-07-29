#include "stdlib.h"
#include <time.h>
#include <stdio.h>
#include "dynamics.h"
#include <math.h>

#define PASO 0.0001
#define PASO2 0.000000005          // esto es paso al cuadrado sobre dos
#define R_CUT2 6.25   // el posta es 6.25

// pasos a seguir:
// calcular la energia (chequear que sea cte)
// calcular presion

int main(int argc, char **argv) {

  FILE *fdat;
  fdat = fopen("../data/md_1a_paso_1e-4.csv", "w");
  fprintf(fdat, "K  U  E\n");
  // Inicializamos numero de particulas
  int n_part = 512;
  float rho = 0.8442; // la posta es 0.8442;
  float l;
  l = pow((float)n_part/rho, 1./3.);
  int steps = 10000;
  float temp_inicial = 0.728;
  srand(time(NULL));

  // Inicializamos arrays de energias cinetica, potencial y total
  float* potential = calloc(steps, sizeof(float));
  float* kinetic = calloc(steps, sizeof(float));
  float* E_total = calloc(steps, sizeof(float));

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

  // inicializamos posiciones, velocidades y fuerzas
  initizalize_pos(n_part, l, pos_x_ant, pos_y_ant, pos_z_ant);
  initizalize_vel(n_part, vel_x_ant, vel_y_ant, vel_z_ant, temp_inicial);
  F_todas(n_part, l, pos_x_ant, pos_y_ant, pos_z_ant, R_CUT2, fuerza_x_ant, fuerza_y_ant, fuerza_z_ant);

  // inicializamos energias cinetica, potencial y total
  kinetic_energy(n_part, 0, vel_x_ant, vel_y_ant, vel_z_ant, kinetic);
  potential_energy(n_part, 0, l, pos_x_ant, pos_y_ant, pos_z_ant, R_CUT2, potential);
  total_energy(0, kinetic, potential, E_total);

  printf("------------------------------------");
  printf("INICIAL\n");
  printf("ENERGIAS\n");
  printf("cinetica = %f, potencial = %f, total = %f\n", kinetic[0], potential[0], E_total[0]);
  printf("\n");
  fprintf(fdat, "%.4g,%.4g,%.4g\n", kinetic[0], potential[0], E_total[0]);

  // iteracion temporal
  int i;
  for (i = 1; i < steps; i++) {

	// se hace la evolucion temporal y se calculan las nuevas energias
    time_evol(n_part, l, PASO, PASO2, pos_x_ant, pos_x_post,
              pos_y_ant, pos_y_post, pos_z_ant,
              pos_z_post, vel_x_ant, vel_x_post,
              vel_y_ant, vel_y_post, vel_z_ant,
              vel_z_post, fuerza_x_ant, fuerza_x_post,
              fuerza_y_ant, fuerza_y_post, fuerza_z_ant,
              fuerza_z_post, R_CUT2);
    kinetic_energy(n_part, i, vel_x_post, vel_y_post, vel_z_post, kinetic);
    potential_energy(n_part, i, l, pos_x_post, pos_y_post, pos_z_post, R_CUT2, potential);
    total_energy(i, kinetic, potential, E_total);


    printf("------------------------------------");
    printf("i = %i\n", i);
    printf("ENERGIAS\n");
    printf("cinetica = %f, potencial = %f, total = %f\n", kinetic[i], potential[i], E_total[i]);
    printf("\n");
    fprintf(fdat, "%.4g,%.4g,%.4g\n", kinetic[i], potential[i], E_total[i]);

    //ahora las post son las nuevas ant para la proxima iteracion
    int j;
    for (j = 0; j < n_part; j++) {
      pos_x_ant[j] = pos_x_post[j];
      pos_y_ant[j] = pos_y_post[j];
      pos_z_ant[j] = pos_z_post[j];

      vel_x_ant[j] = vel_x_post[j];
      vel_y_ant[j] = vel_y_post[j];
      vel_z_ant[j] = vel_z_post[j];

      fuerza_x_ant[j] = fuerza_x_post[j];
      fuerza_y_ant[j] = fuerza_y_post[j];
      fuerza_z_ant[j] = fuerza_z_post[j];
    }

  }

  free(pos_x_ant);
  free(pos_x_post);
  free(pos_y_ant);
  free(pos_y_post);
  free(pos_z_ant);
  free(pos_z_post);
  free(vel_x_ant);
  free(vel_x_post);
  free(vel_y_ant);
  free(vel_y_post);
  free(vel_z_ant);
  free(vel_z_post);
  free(fuerza_x_ant);
  free(fuerza_x_post);
  free(fuerza_y_ant);
  free(fuerza_y_post);
  free(fuerza_z_ant);
  free(fuerza_z_post);
  free(potential);
  free(kinetic);
  free(E_total);


  return 0;
}
