#include "dynamics.h"
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

int hola(){

  printf("hola\n");
  for (int i = 0; i < 0; i++) {
    printf("chau\n");
  }
  return 0;
}


// da valores iniciales a la posicion y velocidad de todas las part
int initizalize_pos(int n_part, float l, float* pos_x_ant, float* pos_y_ant, float* pos_z_ant){
  float a = l/pow(n_part, 1/3.);       // distancia entre particulas en fc del lado de la caja
  // printf("a = %f\n", a);
  int n = 0;
  float b = pow(n_part, 1/3.) - 0.5;
  for (int i = 0; i < b; i++) {
    for (int j = 0; j < b; j++) {
      for (int k = 0; k < b; k++) {
        pos_x_ant[n] = a/2 + i*a;
        // printf("pos_x_ant(%i) = %f\n", n, pos_x_ant[n]);
        pos_y_ant[n] = a/2 + j*a;
        // printf("pos_y_ant(%i) = %f\n", n, pos_y_ant[n]);
        pos_z_ant[n] = a/2 + k*a;
        // printf("pos_y_ant(%i) = %f\n", n, pos_z_ant[n]);
        // printf("\n");
        n++;
      }
    }
  }
  return 0;
}



float gauss(float T){                              // faltaria que dependa de sigma (mu = 0)
  float sigma_gauss2 = 3*T;
  float delta = sqrt(1200*sigma_gauss2);
  printf("delta = %f\n", delta);
  float s = 0;
  for (int i = 0; i < 100; i++) {
    if (i == 1) {
      printf("s[0] = %f\n", s);
    }
    s += ((float) fmodf(rand(), delta));          // habia que agregar (float)
    // printf("paso i =%i, s = %f\n", i, s);
  }
  // printf("gauss() = %f\n", (s/100 -0.5));
  return s/100 -0.5;
}

int initizalize_vel(int n_part, float* vel_x_ant, float* vel_y_ant, float* vel_z_ant, float T){
  for (int n = 0; n < n_part; n++) {
    vel_x_ant[n] = gauss(T);
    vel_y_ant[n] = gauss(T);
    vel_z_ant[n] = gauss(T);
  }
  return 0;
}

// avanza la posicion de una particula dada la fuerza
// r(t+h) = r(t) + hv(t) + h^2/2 F(t)
int adv_pos(int n, float* pos_ant, float* pos_post, float l, float* vel, float paso, float paso2, float* fuerza){

  pos_post[n] = pos_ant[n] + paso*(vel[n]) + (paso2)*(fuerza[n]);
  pos_post[n] = pos_post[n] - l * (fmodf(pos_post[n], l));
  return 0;
}

// avanza la velocidad de una particula con verlet
// v(t+h) = v(t) + h*[F(t+h)+F(t)]/2
int adv_vel(int n, float* vel_ant, float* vel_post, float paso, float* fuerza_ant, float* fuerza_post){

  vel_post[n] = vel_ant[n] + paso*( fuerza_post[n] + fuerza_ant[n] )/2;

  return 0;
}

int sign(float a, float b){
  if (b >= 0) {
    return a;
  }
  else{
    return -a;
  }
}

float delta_coord(int i, int j, float* coord, float l){
  float delta = coord[i] - coord[j];
  if (abs(delta) > 0.5 * l) {
    delta = delta - sign(l, delta);
  }
  return delta;
}

// calcula la distancia entre dos particulas
float dist2(float dx, float dy, float dz){
  float r_ij = ( pow(dx,2) + pow(dy,2) + pow(dz,2) );
  return r_ij + 0.000000001; // para evitar que sea cero, ver el orden de la correccion
}

// calcula el factor de las fuerza (usando derivada de LJ) entre dos particulas
// separadas por una distancia/ (distancia)^2 = dist2
float eval_LJ(float dist2, float r_cut2){      // sigma como macro? unidades reducidas?
  // float epsilon = 1;
  // float sigma = 1;
  // float sigma6 = pow(sigma, 6);
  float invr2 = 1/dist2;
  float F = 0;
  if (dist2 < r_cut2) {
    // F = 24 * epsilon * sigma6 * pow(invr2, 3) * (2 * sigma6 * pow(invr2, 4) - invr2);
    F = 48 * pow(invr2, 4) * (pow(invr2, 3) - 0.5);
  }
  return F;
}


int F_tot(int i, int n_part, float l, float* x, float* y, float* z, float r_cut2, float *fuerza_x, float *fuerza_y, float *fuerza_z){
  // vacia el array de fuerzas antes de calcular
  fuerza_x[i] = 0;
  fuerza_y[i] = 0;
  fuerza_z[i] = 0;
  float distancia2 = 0;
  float F = 0;
  int j;
  // suma para cada particula i, las fuerzas F_ij, que resulta en la F_i total

  float dx, dy, dz;
  for (j = 0; j < i; j++) {
    dx = delta_coord(i, j, x, l);
    dy = delta_coord(i, j, y, l);
    dz = delta_coord(i, j, z, l);
    distancia2 = dist2(dx, dy, dz);
    F = eval_LJ(distancia2, r_cut2);
    // printf("F de particula i = %i ejercida por j = %i, distancia2 =  %f, \n coord x_i = %f, coord y_i = %f, coord z_i = %f, \n coord x_j = %f, coord y_j = %f, coord z_j = %f, \n F = %f\n", i, j, distancia2, x[i], y[i], z[i], x[j], y[j], z[j], F);
    fuerza_x[i] += F * (dx);
    fuerza_y[i] += F * (dy);
    fuerza_z[i] += F * (dz);
  }
  // lo divide en rangos para evitar i=j
  for (j = i+1; j < n_part; j++) {
    dx = delta_coord(i, j, x, l);
    dy = delta_coord(i, j, y, l);
    dz = delta_coord(i, j, z, l);
    distancia2 = dist2(dx, dy, dz);
    F = eval_LJ(distancia2, r_cut2);
    // printf("(fuerza sobre la particula i = %i ejercida por j = %i con una distancia2 = %f) F = %f\n", i, j, distancia2, F);
    fuerza_x[i] += F * (dx);
    fuerza_y[i] += F * (dy);
    fuerza_z[i] += F * (dz);
  }
  return 0;
}

int time_evol(int n_part, float l, float paso, float paso2, float* pos_x_ant, float* pos_x_post,
              float* pos_y_ant, float* pos_y_post, float* pos_z_ant,
              float* pos_z_post, float* vel_x_ant, float* vel_x_post,
              float* vel_y_ant, float* vel_y_post, float* vel_z_ant,
              float* vel_z_post, float* fuerza_x_ant, float* fuerza_x_post,
              float* fuerza_y_ant, float* fuerza_y_post, float* fuerza_z_ant,
              float* fuerza_z_post, float r_cut2){
  for (int n = 0; n < n_part; n++) {
    adv_pos(n, pos_x_ant, pos_x_post, l, vel_x_ant, paso, paso2, fuerza_x_ant);
    adv_pos(n, pos_y_ant, pos_y_post, l, vel_y_ant, paso, paso2, fuerza_y_ant);
    adv_pos(n, pos_z_ant, pos_z_post, l, vel_z_ant, paso, paso2, fuerza_z_ant);

    F_tot(n, n_part, l, pos_x_post, pos_y_post, pos_z_post, r_cut2, fuerza_x_post, fuerza_y_post, fuerza_z_post);

    adv_vel(n, vel_x_ant, vel_x_post, paso, fuerza_x_ant, fuerza_x_post);
    adv_vel(n, vel_y_ant, vel_y_post, paso, fuerza_y_ant, fuerza_y_post);
    adv_vel(n, vel_z_ant, vel_z_post, paso, fuerza_z_ant, fuerza_z_post);
  }
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
