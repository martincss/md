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
  fdat = fopen("../data/md_3_prueba_posta_3.csv", "w");
  fprintf(fdat, "r  g1  g2  g3\n");
  // Inicializamos numero de particulas
  int n_part = 512;
  float rho = 0.6; // la posta es 0.8442;
  float l;
  l = pow((float)n_part/rho, 1./3.);
  int steps_term = 4000;
  int tiempo_desc = 500;
  float temp_inicial = 1.5;
  srand(time(NULL));

  float delta_r = 0.025;
  int particiones_r = (int)(l*0.5/delta_r);
  int cant_sample = 10;

  // inicializamos array de r, y de las g
  float* radios = calloc(particiones_r, sizeof(float)); // tiene el tamano de la cant de r
  float* g_sample = calloc(particiones_r*cant_sample, sizeof(float)); // tiene la cant de r, tantas veces como la cant_sampleos
  float* g_prom = calloc(particiones_r*3, sizeof(float)); // tiene la cant de r, tantas como valores de rho (tres)

  for (int k = 0; k < particiones_r; k++) {
    radios[k] = k*delta_r;
  }

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

  rho = 0.4;
  for (int r = 0; r < 3; r++) {
    rho += 0.2;
	l = pow((float)n_part/rho, 1./3.);

	// vacía el g_sample
	for(int a = 0; a < particiones_r*cant_sample; a++){
		g_sample[a] = 0;
	}
	
    printf("========================================== RHO = %f\n", rho);
    // inicializamos posiciones, velocidades y fuerzas
    initizalize_pos(n_part, l, pos_x_ant, pos_y_ant, pos_z_ant);
    initizalize_vel(n_part, vel_x_ant, vel_y_ant, vel_z_ant, temp_inicial);
    F_todas(n_part, l, pos_x_ant, pos_y_ant, pos_z_ant, R_CUT2, fuerza_x_ant, fuerza_y_ant, fuerza_z_ant);

    // iteracion temporal hasta termalizar
    int i;
    for (i = 0; i < steps_term; i++) {

  	// se hace la evolucion temporal y se calculan las nuevas energias
      time_evol(n_part, l, PASO, PASO2, pos_x_ant, pos_x_post,
                pos_y_ant, pos_y_post, pos_z_ant,
                pos_z_post, vel_x_ant, vel_x_post,
                vel_y_ant, vel_y_post, vel_z_ant,
                vel_z_post, fuerza_x_ant, fuerza_x_post,
                fuerza_y_ant, fuerza_y_post, fuerza_z_ant,
                fuerza_z_post, R_CUT2);

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

    // hace el sampleo
    for (i = 0; i < cant_sample*tiempo_desc; i++) {

  	// se hace la evolucion temporal y se calculan las nuevas energias
      time_evol(n_part, l, PASO, PASO2, pos_x_ant, pos_x_post,
                pos_y_ant, pos_y_post, pos_z_ant,
                pos_z_post, vel_x_ant, vel_x_post,
                vel_y_ant, vel_y_post, vel_z_ant,
                vel_z_post, fuerza_x_ant, fuerza_x_post,
                fuerza_y_ant, fuerza_y_post, fuerza_z_ant,
                fuerza_z_post, R_CUT2);

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

      // cuando el tiempo es multiplo del tiempo de sampleo, calcula la g para todos los r
      if (i % tiempo_desc == 0) {
        printf("dato tomado para %i de %i con rho = %f\n", i/(tiempo_desc)+1, cant_sample, rho);
        // for (int k = 0; k < particiones_r; k++) {
          // g_sample[ (i/tiempo_desc)*particiones_r + k ] = g(rho, l, radios[k], delta_r, n_part, pos_x_post, pos_y_post, pos_z_post);
        // }
		g_posta(delta_r, l, rho, n_part, radios, g_sample, (i/tiempo_desc)*particiones_r, pos_x_post, pos_y_post, pos_z_post);
      }

    }
    // calcula la g(r) promedio para cada rho
    for (int s = 0; s < cant_sample; s++) { // suma sobre los sampleos
      for (int k = 0; k < particiones_r; k++) { // calcula para cada r
        g_prom[k + r*particiones_r] += g_sample[ s*particiones_r + k]/cant_sample;
      }
    }

  }

  // imprime todos los resultados
  printf("================================ RESULTADOS FINALES\n");
  for (int k = 0; k < particiones_r; k++) { // printea todas para cada r
    printf("r = %f g_1 = %f, g_2 = %f, g_3 = %f\n", radios[k], g_prom[k], g_prom[k + particiones_r], g_prom[k + 2*particiones_r]);
    fprintf(fdat, "%.3g,%.6g,%.6g,%.6g\n", radios[k], g_prom[k], g_prom[k + particiones_r], g_prom[k + 2*particiones_r]);
  }
  fflush(fdat);
  fclose(fdat);

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


  return 0;
}
