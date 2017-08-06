#include "stdlib.h"
#include <time.h>
#include <stdio.h>
#include "dynamics.h"
#include <math.h>

#define PASO 0.0001
#define PASO2 0.000000005          // esto es paso al cuadrado sobre dos
#define R_CUT2 6.25   // el posta es 6.25

int main(int argc, char **argv) {

  FILE *fdat;
  fdat = fopen("../data/md_2_prueba.csv", "w");
  fprintf(fdat, "rho  temp  E  P\n");

  // Inicializamos numero de particulas
  srand(time(NULL));
  int n_part = 125;
  float l;
  int steps_term = 4000; // tiempo de termalizacion
  int tiempo_desc = 500; // tiempo de descorrelacion
  //=======================================================
  float rho_inicial = 0.4; // tomar acá el valor mínimo, para que empieze por la l máxima
  float rho_final = 0.8;
  float temp_inicial = 2.0;
  float temp_final = 0.4;

  int cant_rho = 5;							// cantidad de densidades
  int cant_temps = 10;           // cantidad de temperaturas
  int cant_sample = 10;					// cantidad de samples MENOS UNO (por como se da el loop)

  float incr_temp = (temp_inicial - temp_final)/(cant_temps-1);
  float incr_rho = (rho_final - rho_inicial)/(cant_rho-1);

  // Inicializamos pointers de sampleo
  float* potential = calloc(1, sizeof(float));
  float* kinetic = calloc(1, sizeof(float));
  float* E_sample = calloc(1, sizeof(float));
  float P_sample;
  float temperatura = temp_inicial;
  float rho = rho_inicial;
  float temperatura_samp;

  // Inicializamos pointers de resultados promedio
  float* E_prom = calloc(cant_rho*cant_temps, sizeof(float));
  float* p_prom = calloc(cant_rho*cant_temps, sizeof(float));


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

  // INICIA EL LOOP EN RHO
  for (int r = 0; r < cant_rho; r++) {

    l = pow((float)n_part/rho, 1./3.); // debe redefinirse para cada rho


    printf("========================================== RHO = %f\n", rho);
    printf("================================= TEMP = %f\n", temp_inicial);
    /*
    Tomamos el sample para la temperatura inicial y luego hacemos el
    loop sobre temperaturas con rescaling, tomando datos.
    */
    // inicializamos posiciones, velocidades y fuerzas
    initizalize_pos(n_part, l, pos_x_ant, pos_y_ant, pos_z_ant);
    initizalize_vel(n_part, vel_x_ant, vel_y_ant, vel_z_ant, temp_inicial);
    F_todas(n_part, l, pos_x_ant, pos_y_ant, pos_z_ant, R_CUT2, fuerza_x_ant, fuerza_y_ant, fuerza_z_ant);

    // iteracion temporal hasta termalizar
	// =========================================================
    int i;
    for (i = 0; i < steps_term; i++) {

      time_evol(n_part, l, PASO, PASO2, pos_x_ant, pos_x_post,
                pos_y_ant, pos_y_post, pos_z_ant,
                pos_z_post, vel_x_ant, vel_x_post,
                vel_y_ant, vel_y_post, vel_z_ant,
                vel_z_post, fuerza_x_ant, fuerza_x_post,
                fuerza_y_ant, fuerza_y_post, fuerza_z_ant,
                fuerza_z_post, R_CUT2);

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
    printf("\n HECHA LA TERMALIZACION\n \n");
	// ========================================================

	/* a partir de acá, se evoluciona el sistema desde termalizacion,
	y cada tiempo_desc pasos se toman datos para el sampleo */
    for (i = 0; i < cant_sample*tiempo_desc + 1; i++) {
	  // evolucion temporal ====================================
      time_evol(n_part, l, PASO, PASO2, pos_x_ant, pos_x_post,
                pos_y_ant, pos_y_post, pos_z_ant,
                pos_z_post, vel_x_ant, vel_x_post,
                vel_y_ant, vel_y_post, vel_z_ant,
                vel_z_post, fuerza_x_ant, fuerza_x_post,
                fuerza_y_ant, fuerza_y_post, fuerza_z_ant,
                fuerza_z_post, R_CUT2);

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
      // printf("iteracion post-termalizacion i = %i\n", i);
	  // =======================================================

      // cuando el tiempo es multiplo del tiempo de sampleo, calcula la g para todos los r
      if (i % tiempo_desc == 0) {
        printf("dato tomado para %i de %i con rho = %f, T = %f\n", i/(tiempo_desc)+1, cant_sample, rho, temp_inicial);
        kinetic_energy(n_part, 0, vel_x_post, vel_y_post, vel_z_post, kinetic);
        potential_energy(n_part, 0, l, pos_x_post, pos_y_post, pos_z_post, R_CUT2, potential);
        total_energy(0, kinetic, potential, E_sample);
        temperatura_samp = temperature(n_part, vel_x_post, vel_y_post, vel_z_post);
        P_sample = pressure(n_part, l, rho, pos_x_post, pos_y_post, pos_z_post, R_CUT2, temperatura_samp);
        printf("con t_samp = %f tenemos presion = %f\n", temperatura_samp, P_sample);
        E_prom[r*cant_temps] += E_sample[0]/(cant_sample+1);
        p_prom[r*cant_temps] += P_sample/(cant_sample+1);
      }
    }

    // EMPIEZA EL LOOP EN TEMPERATURAS, DENTRO DEL LOOP DE RHOS (o sea, para cada rho)
    temperatura = temp_inicial;
    for (int s = 1; s < cant_temps; s++) {
      temperatura -= incr_temp;
      printf("===================================== TEMP = %f\n", temperatura);
      rescaling(temperatura_samp, temperatura, n_part, vel_x_ant, vel_y_ant, vel_z_ant);
      printf("HECHO EL RESCALING \n");
      // iteracion temporal hasta termalizar
    // =========================================================
      int i;
      for (i = 0; i < steps_term; i++) {

        time_evol(n_part, l, PASO, PASO2, pos_x_ant, pos_x_post,
                  pos_y_ant, pos_y_post, pos_z_ant,
                  pos_z_post, vel_x_ant, vel_x_post,
                  vel_y_ant, vel_y_post, vel_z_ant,
                  vel_z_post, fuerza_x_ant, fuerza_x_post,
                  fuerza_y_ant, fuerza_y_post, fuerza_z_ant,
                  fuerza_z_post, R_CUT2);

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
      printf("\n HECHA LA TERMALIZACION\n \n");
    // ========================================================

    /* a partir de acá, se evoluciona el sistema desde termalizacion,
    y cada tiempo_desc pasos se toman datos para el sampleo */
      for (i = 0; i < cant_sample*tiempo_desc + 1; i++) {
      // evolucion temporal ====================================
        time_evol(n_part, l, PASO, PASO2, pos_x_ant, pos_x_post,
                  pos_y_ant, pos_y_post, pos_z_ant,
                  pos_z_post, vel_x_ant, vel_x_post,
                  vel_y_ant, vel_y_post, vel_z_ant,
                  vel_z_post, fuerza_x_ant, fuerza_x_post,
                  fuerza_y_ant, fuerza_y_post, fuerza_z_ant,
                  fuerza_z_post, R_CUT2);

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
      // =======================================================

        // cuando el tiempo es multiplo del tiempo de sampleo, calcula la g para todos los r
        if (i % tiempo_desc == 0) {
          printf("dato tomado para %i de %i con rho = %f, T = %f\n", i/(tiempo_desc)+1, cant_sample, rho, temperatura);
          kinetic_energy(n_part, 0, vel_x_post, vel_y_post, vel_z_post, kinetic);
          potential_energy(n_part, 0, l, pos_x_post, pos_y_post, pos_z_post, R_CUT2, potential);
          total_energy(0, kinetic, potential, E_sample);
          temperatura_samp = temperature(n_part, vel_x_post, vel_y_post, vel_z_post);
          P_sample = pressure(n_part, l, rho, pos_x_post, pos_y_post, pos_z_post, R_CUT2, temperatura_samp);
          printf("con t_samp = %f tenemos presion = %f\n", temperatura_samp, P_sample);
          // // s es la iteracion actual de temperatura
          E_prom[r*cant_temps + s] += E_sample[0]/(cant_sample + 1);
          p_prom[r*cant_temps + s] += P_sample/(cant_sample + 1);
        }
      }

    // ACA TERMINA EL LOOP EN TEMPERATURAS (para cada rho)
    }
    rho += incr_rho; // aumenta el rho acá para que el primero sea el incial
  // ACA TERMINA EL LOOP EN RHOS
  }

  // imprime todos los resultados y los guarda al archivo
  printf("================================ RESULTADOS FINALES\n");

  rho = rho_inicial - incr_rho;
  temperatura = temp_inicial + incr_temp;
  for (int r = 0; r < cant_rho; r++) {
      rho += incr_rho;
      temperatura = temp_inicial + incr_temp;
      for (int s = 0; s < cant_temps; s++) {
        temperatura -= incr_temp;
        printf("rho = %.4f temp = %.4f E = %.4f P = %.4f\n", rho, temperatura, E_prom[r*cant_temps + s], p_prom[r*cant_temps + s]);
        fprintf(fdat,"%.4g,%.4g,%.4g,%.4g\n", rho, temperatura, E_prom[r*cant_temps + s], p_prom[r*cant_temps + s]);
      }
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
  free(potential);
  free(kinetic);
  free(E_sample);
  free(E_prom);
  free(p_prom);

  return 0;
}
