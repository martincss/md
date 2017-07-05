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
  int n_part = 8;
  float l = 2;
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

//-----------------------------------------------------------------------------

  initizalize_pos(n_part, l, pos_x_ant, pos_y_ant, pos_z_ant);
  initizalize_vel(n_part, vel_x_ant, vel_y_ant, vel_z_ant, T);

  kinetic_temperature(n_part, 0, vel_x_ant, vel_y_ant, vel_z_ant, temperature, kinetic);
  // potential_energy(n_part, 0, l, pos_x_ant, pos_y_ant, pos_z_ant, R_CUT2, potential);
  // total_energy(0, kinetic, potential, E_total);

  printf("------------------------------------");
  printf("INICIAL\n");
  // printf("ENERGIAS\n");
  // printf("cinetica = %f, potencial = %f, total = %f\n", kinetic[0], potential[0], E_total[0]);
  printf("temperatura = %f\n", temperature[0]);
  printf("\n");

//-----------------------------------------------------------------------------

  int i;
  for (i = 0; i < steps; i++) {

  // int k;
  // for (k = 0; k < n_part; k++) {
  //   printf("particula = %i,   ", k);
  //   printf("\n");
  //   printf("pos x es %f,    ", pos_x_ant[k]);
  //   printf("vel x es %f,    ", vel_x_ant[k]);
  //   printf("fuerza x es %f,    ", fuerza_x_ant[k]);
  //   printf("\n");
  //   printf("pos y es %f,    ", pos_y_ant[k]);
  //   printf("vel y es %f,    ", vel_y_ant[k]);
  //   printf("fuerza y es %f,    ", fuerza_y_ant[k]);
  //   printf("\n");
  //   printf("pos z es %f,    ", pos_z_ant[k]);
  //   printf("vel z es %f,    ", vel_z_ant[k]);
  //   printf("fuerza z es %f,    ", fuerza_z_ant[k]);
  //   printf("\n");
  //   printf("\n");
  // }
  // esto es el time
  for (int n = 0; n < n_part; n++) {
    adv_pos(n, pos_x_ant, pos_x_post, l, vel_x_ant, PASO, PASO2, fuerza_x_ant);
    adv_pos(n, pos_y_ant, pos_y_post, l, vel_y_ant, PASO, PASO2, fuerza_y_ant);
    adv_pos(n, pos_z_ant, pos_z_post, l, vel_z_ant, PASO, PASO2, fuerza_z_ant);
    printf("particula %i, pos_x_ant = %f\n", n, pos_x_ant[n]);
    printf("particula %i, pos_y_ant = %f\n", n, pos_y_ant[n]);
    printf("particula %i, pos_z_ant = %f\n", n, pos_z_ant[n]);
    printf("\n");
    printf("particula %i, pos_x_post = %f\n", n, pos_x_post[n]);
    printf("particula %i, pos_y_post = %f\n", n, pos_y_post[n]);
    printf("particula %i, pos_z_post = %f\n", n, pos_z_post[n]);


    F_tot(n, n_part, l, pos_x_post, pos_y_post, pos_z_post, R_CUT2, fuerza_x_post, fuerza_y_post, fuerza_z_post);

    // adv_vel(n, vel_x_ant, vel_x_post, PASO, fuerza_x_ant, fuerza_x_post);
    // adv_vel(n, vel_y_ant, vel_y_post, PASO, fuerza_y_ant, fuerza_y_post);
    // adv_vel(n, vel_z_ant, vel_z_post, PASO, fuerza_z_ant, fuerza_z_post);
  }

    // time_evol(n_part, l, PASO, PASO2, pos_x_ant, pos_x_post,
    //           pos_y_ant, pos_y_post, pos_z_ant,
    //           pos_z_post, vel_x_ant, vel_x_post,
    //           vel_y_ant, vel_y_post, vel_z_ant,
    //           vel_z_post, fuerza_x_ant, fuerza_x_post,
    //           fuerza_y_ant, fuerza_y_post, fuerza_z_ant,
    //           fuerza_z_post, R_CUT2);
    printf("---------------------------------time evoluciono \n");
    // kinetic_temperature(n_part, i, vel_x_post, vel_y_post, vel_z_post, temperature, kinetic);
    // potential_energy(n_part, i, l, pos_x_post, pos_y_post, pos_z_post, R_CUT2, potential);
    // total_energy(i, kinetic, potential, E_total);
    // printf("i = %i\n", i);
    // printf("ENERGIAS\n");
    // printf("cinetica = %f, potencial = %f, total = %f\n", kinetic[i], potential[i], E_total[i]);
    // printf("temperatura = %f\n", temperature[i]);
    // printf("\n");

    //ahora las post son las nuevas ant
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

    printf("-------------------------------PRUEBA DIST2 \n");
    prueba_dist(n_part, l, pos_x_ant, pos_y_ant, pos_z_ant);

  }
// printf("------------------------------------ F TOT AISLADA");
// for (i = 0; i < n_part; i++) {
//   F_tot(i, n_part, l, pos_x_ant, pos_y_ant, pos_z_ant, R_CUT2, fuerza_x_post, fuerza_y_post, fuerza_z_post);
// }

  free(pos_x_ant);
  free(pos_x_post);
  free(pos_y_ant);
  free(pos_y_post);

  return 0;
}
