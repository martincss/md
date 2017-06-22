#include "dynamics.h"
#include <time.h>
#include <stdlib.h>
#include <math.h>

int hola(){

  printf("hola\n");

  return 0;
}

// da valores iniciales a la posicion y velocidad de todas las part
int initizalize(int n_part, float* pos_x_ant, float* pos_y_ant, float* pos_z_ant,
                float* vel_x_ant, float* vel_y_ant, float* vel_z_ant){
    // *** le puse cero por ahora, despues ver
    for (int i = 0; i < n_part; i++) {
      pos_x_ant[i] = 0;
      pos_y_ant[i] = 0;
      pos_z_ant[i] = 0;
      vel_x_ant[i] = 0;
      vel_y_ant[i] = 0;
      vel_z_ant[i] = 0;
    }

}

// avanza la posicion de una particula dada la fuerza
// r(t+h) = r(t) + hv(t) + h^2/2 F(t)
int adv_pos(float* pos_ant, float* pos_post, float* vel, float paso, float* fuerza){

  *pos_post = *pos_ant + paso*(*vel) + (paso*paso/2)*(*fuerza);

  return 0;
}

// avanza la velocidad de una particula con verlet
// v(t+h) = v(t) + h*[F(t+h)+F(t)]/2
int adv_vel(float* vel_ant, float* vel_post, float paso, float* fuerza_ant, float* fuerza_post){

  *vel_post = *vel_ant + paso*( *fuerza_post + *fuerza_ant )/2;

  return 0;
}

// calcula la temperatura a partir del valor medio de la energia cinetica
float temperature(int n_part, float* vel_x, float* vel_y, float* vel_z){

    float KE = 0, temp = 0;
    // calcula la energia cinetica total (KE)
    for (size_t i = 0; i < n_part; i++) {
        KE += (vel_x[i]*vel_x[i] + vel_y[i]*vel_y[i] + vel_z[i]*vel_z[i]);
    }
    KE = KE/2;
    temp = (2/3)*KE/n_part;

    return temp;
}
