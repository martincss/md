#include "stdlib.h"
#include <time.h>
#include <stdio.h>
#include "dynamics.h"
#include "pruebas_pos.h"
#include <math.h>

#define PASO 0.01
#define PASO2 0.00005          // esto es paso al cuadrado sobre dos
#define R_CUT2 6.25   // el posta es 6.25

int main(int argc, char **argv) {

  hola();
  // Inicializamos numero de particulas
  int n_part = 216;
  float rho = 0.4; // la posta es 0.8442;
  float l;
  l = pow((float)n_part/rho, 1./3.);
  int steps = 100;
  int t_samp = 1;
  float temp_inicial = 10.0;
  srand(time(NULL));

  /*
  Abrimos el archivo de registro de posiciones
  En tres líneas consecutivas se registrarán las posiciones x y z de
  todas las partículas, para un instante temporal.
  */
  FILE *fdat;
  fdat = fopen("../data/animation/n216_T100_t100.csv", "w");
  fprintf(fdat, "n = %i\n", n_part);

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

  // init_prueba_dos(l, pos_x_ant, pos_y_ant, pos_z_ant);
  // init_prueba_tres(l, pos_x_ant, pos_y_ant, pos_z_ant);
  // init_prueba_cuatro(l, pos_x_ant, pos_y_ant, pos_z_ant);
  // init_prueba_cinco(l, pos_x_ant, pos_y_ant, pos_z_ant);

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
  int k;
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
  printf("ENERGIAS\n");
  printf("cinetica = %f, potencial = %f, total = %f\n", kinetic[0], potential[0], E_total[0]);
  printf("\n");
  // sum_vel(n_part, vel_x_ant, vel_y_ant, vel_z_ant);


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

    // for (int n = 0; n < n_part; n++) {
    //   printf("la fuerza_x sobre n=%i es %f \n", n, fuerza_x_post[n]);
    //   printf("la fuerza_y sobre n=%i es %f \n", n, fuerza_y_post[n]);
    //   printf("la fuerza_z sobre n=%i es %f \n", n, fuerza_z_post[n]);
    // }

    printf("------------------------------------");
    printf("i = %i\n", i);
  //
	// int k;
  // for (k = 0; k < n_part; k++) {
  //   printf("particula = %i,   ", k);
  //   printf("\n");
  //   printf("pos x es %f,    ", pos_x_post[k]);
  //   printf("vel x es %f,    ", vel_x_post[k]);
  //   printf("fuerza x es %f,    ", fuerza_x_post[k]);
  //   printf("\n");
  //   printf("pos y es %f,    ", pos_y_post[k]);
  //   printf("vel y es %f,    ", vel_y_post[k]);
  //   printf("fuerza y es %f,    ", fuerza_y_post[k]);
  //   printf("\n");
  //   printf("pos z es %f,    ", pos_z_post[k]);
  //   printf("vel z es %f,    ", vel_z_post[k]);
  //   printf("fuerza z es %f,    ", fuerza_z_post[k]);
  //   printf("\n");
  //   printf("\n");
  // }

    printf("ENERGIAS\n");
    printf("cinetica = %f, potencial = %f, total = %f\n", kinetic[i], potential[i], E_total[i]);
    printf("\n");
    // printf("%f, %f, %f\n", kinetic[i], potential[i], kinetic[i]+potential[i]);
	// sum_vel(n_part, vel_x_post, vel_y_post, vel_z_post);

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

	if (i % t_samp == 0){
		for(j = 0; j < n_part-1; j++){
			fprintf(fdat, "%.3g,", pos_x_post[j]);
		}
		fprintf(fdat, "%.3g", pos_x_post[n_part-1]);
		fprintf(fdat, "\n");
		for(j = 0; j < n_part-1; j++){
			fprintf(fdat, "%.3g,", pos_y_post[j]);
		}
		fprintf(fdat, "%.3g", pos_y_post[n_part-1]);
		fprintf(fdat, "\n");
		for(j = 0; j < n_part-1; j++){
			fprintf(fdat, "%.3g,", pos_z_post[j]);
		}
		fprintf(fdat, "%.3g", pos_z_post[n_part-1]);
		fprintf(fdat, "\n");
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
