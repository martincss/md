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



// da valores iniciales a la posicion de todas las part ubicadas en un retículo fcc
int initizalize_pos(int n_part, float l, float* pos_x_ant, float* pos_y_ant, float* pos_z_ant){
  float a = l/pow(n_part, 1/3.);       // distancia entre particulas en fc del lado de la caja
  // printf("a = %f\n", a);
  int n = 0;
  float b = pow(n_part, 1/3.) - 1;
  for (int i = 0; i <= b; i++) {
    for (int j = 0; j <= b; j++) {
      for (int k = 0; k <= b; k++) {
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

// genera numeros aleatorios con dist normal, con dispersión según la temperatura deseada
float gauss(){                              // faltaria que dependa de sigma (mu = 0)
  // float sigma_gauss2 = T;
  // float delta = sqrt(1200*sigma_gauss2);
  // printf("delta = %f\n", delta);
  float s = 0;
  for (int i = 0; i < 100; i++) {
    s += ((float) rand() / RAND_MAX);          // habia que agregar (float)
  }
  // printf("gauss() = %f\n", (s/100 -0.5));
  // printf("gauss = %f\n", (s/100 -0.5) * delta);
  return (s/100 -0.5);
}

// inicializa velocidades siguiendo una dist. gaussiana para cada componente, y les resta
// el valor medio, para evitar el movimiento neto
int initizalize_vel(int n_part, float* vel_x_ant, float* vel_y_ant, float* vel_z_ant, float T){
  for (int n = 0; n < n_part; n++) {
    vel_x_ant[n] = gauss(T);
    vel_y_ant[n] = gauss(T);
    vel_z_ant[n] = gauss(T);
  }
  // les saca el valor medio..
  float s_x = 0;
  float s_y = 0;
  float s_z = 0;
  for (int n = 0; n < n_part; n++){
	  s_x += vel_x_ant[n];
	  s_y += vel_y_ant[n];
	  s_z += vel_z_ant[n];
  }
  for (int n = 0; n < n_part; n++){
	  vel_x_ant[n] -= s_x/n_part;
	  vel_y_ant[n] -= s_y/n_part;
	  vel_z_ant[n] -= s_z/n_part;
  }

  // calcula la temperatura
  float T_act = 0;
  for (int n = 0; n < n_part; n++) {
    T_act += vel_x_ant[n]*vel_x_ant[n] + vel_y_ant[n]*vel_y_ant[n] + vel_z_ant[n]*vel_z_ant[n];
  }
  T_act /= 3*n_part;

  // rescaling para lograr la deseada
  for (int n = 0; n < n_part; n++){
	  vel_x_ant[n] *= sqrt(T/T_act);
	  vel_y_ant[n] *= sqrt(T/T_act);
	  vel_z_ant[n] *= sqrt(T/T_act);
  }


  return 0;
}
// para verificar que no haya velocidad media no nula
int sum_vel(int n_part, float* vel_x, float* vel_y, float* vel_z){
	float sum_x = 0;
	float sum_y = 0;
	float sum_z = 0;
	int i;
	for (i = 0; i < n_part; i++){
		sum_x += vel_x[i];
		sum_y += vel_y[i];
		sum_z += vel_z[i];
	}
	printf("las velocidades promedio son v_x = %f, v_y = %f, v_z = %f \n \n", sum_x/n_part, sum_y/n_part, sum_z/n_part);
}

// avanza la posicion de una particula dada la fuerza
// r(t+h) = r(t) + hv(t) + h^2/2 F(t)
int adv_pos(int n, float* pos_ant, float* pos_post, float l, float* vel, float paso, float paso2, float* fuerza){

  pos_post[n] = pos_ant[n] + paso*(vel[n]) + (paso2)*(fuerza[n]);
  pos_post[n] = fmodf(pos_post[n],l);
  if (pos_post[n] < 0) {
    pos_post[n] += l;
  }
  if (pos_post[n] > l || pos_post[n] < 0) {
    printf("la posicion se fue del rango\n");
  }
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

// hace la delta de una coordenada entre la part. i y la imagen más cercana de la j (o la j misma)
float delta_coord(int i, int j, float* coord, float l){
  float delta = coord[i] - coord[j];
  while (delta > l/2) delta -= l;
  while (delta < -l/2) delta += l;
  // if (abs(delta) > 0.5 * l) {
  //   delta = delta - sign(l, delta);
  // }
  return delta;
}

// calcula la distancia (al cuadrado) entre dos particulas dadas las deltas de cada coord
float dist2(float dx, float dy, float dz){
  float r_ij = ( pow(dx,2) + pow(dy,2) + pow(dz,2) );
  if (r_ij == 0) {
    printf("se evito que la distancia fuese nula\n");
  }
  if (r_ij < 1){
    ;//printf("========================================== DIST MAS CHICA QUE UNO\n");
  }
  return r_ij;// + 0.000001; // para evitar que sea cero, ver el orden de la correccion
}

// calcula el factor de las fuerza (usando derivada de LJ) entre dos particulas
// separadas por una distancia/ (distancia)^2 = dist2
// le agrego el shift de fuerza para sualizar el cut_off, ver Haile ss. 5.1.2
float eval_LJ(float dist2, float r_cut){

  float invr2 = 1./dist2;
  float F = 0;
  if (dist2 < r_cut) {
    F = 24 * pow(invr2, 3) * (2 * pow(invr2, 4) - invr2);
	//F -= 24 * pow(r_cut, -3) * (2 * pow(r_cut, -4) - 1/r_cut);
	// F = 24 * pow(invr2, 3) * (2 * pow(invr2, 4) - invr2) + 0.015599; // con el shift que es -F(r_cut)

  }
  if (F>10000){
    printf("la fuerza se paso de 10000\n");
  }
  return F;
}

// calcula todas las fuerzas sobre la particula i, (ya no la usamos)
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
  // for (j = 0; j < i; j++) {
  //   dx = delta_coord(i, j, x, l);
  //   dy = delta_coord(i, j, y, l);
  //   dz = delta_coord(i, j, z, l);
  //   distancia2 = dist2(dx, dy, dz);
  //   F = eval_LJ(distancia2, r_cut2);
  //   // printf("F = %f\n", F);
  //   // printf("fuerza ejercida sobre i = %i por j = %i es fuerza = %f \n", i, j, F);
  //   fuerza_x[i] += F * (dx);
  //   fuerza_y[i] += F * (dy);
  //   fuerza_z[i] += F * (dz);
  // }
  // lo divide en rangos para evitar i=j
  for (j = i+1; j < n_part; j++) {
    dx = delta_coord(i, j, x, l);
    dy = delta_coord(i, j, y, l);
    dz = delta_coord(i, j, z, l);
    distancia2 = dist2(dx, dy, dz);
    F = eval_LJ(distancia2, r_cut2);
    // printf("F = %f\n", F);
    // printf("fuerza ejercida sobre i = %i por j = %i es fuerza = %f \n", i, j, F);
    fuerza_x[i] += F * (dx);
    fuerza_y[i] += F * (dy);
    fuerza_z[i] += F * (dz);

    fuerza_x[j] -= F * (dx);
    fuerza_y[j] -= F * (dy);
    fuerza_z[j] -= F * (dz);
  }
  return 0;
}

// calcula todas las fuerzas sobre todas las partículas a la vez
// esta lo hace segun el libro, todas a la vez, no hace falta iterarla
int F_todas(int n_part, float l, float* x, float* y, float* z, float r_cut2, float *fuerza_x, float *fuerza_y, float *fuerza_z){

  float distancia2 = 0;
  float F = 0;
  int i, j;

  float dx, dy, dz;
  for (i = 0; i < n_part-1; i++){
    for (j = i+1; j < n_part; j++) {
      dx = delta_coord(i, j, x, l);
      dy = delta_coord(i, j, y, l);
      dz = delta_coord(i, j, z, l);
      distancia2 = dist2(dx, dy, dz);
      F = eval_LJ(distancia2, r_cut2);
      // printf("fuerza del par i=%i, j=%i\n", i, j);
      // printf("F = %f\n", F);
      // printf("fuerza ejercida sobre i = %i por j = %i es fuerza = %f \n", i, j, F);
      fuerza_x[i] += F * (dx);
      fuerza_y[i] += F * (dy);
      fuerza_z[i] += F * (dz);

      fuerza_x[j] -= F * (dx);
      fuerza_y[j] -= F * (dy);
      fuerza_z[j] -= F * (dz);
    }
  }
  return 0;
}

// genera la evolución temporal de un paso, actualizando secunecialmente posiciones, fuerzas y velocidades
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
  }
  for (int n = 0; n < n_part; n++) {
    fuerza_x_post[n] = 0;
    fuerza_y_post[n] = 0;
    fuerza_z_post[n] = 0;
  }
  F_todas(n_part, l, pos_x_post, pos_y_post, pos_z_post, r_cut2, fuerza_x_post, fuerza_y_post, fuerza_z_post);

  // for (int n = 0; n < n_part-1; n++) {
  //   F_tot(n, n_part, l, pos_x_post, pos_y_post, pos_z_post, r_cut2, fuerza_x_post, fuerza_y_post, fuerza_z_post);
  // }
  for (int n = 0; n < n_part; n++) {
    adv_vel(n, vel_x_ant, vel_x_post, paso, fuerza_x_ant, fuerza_x_post);
    adv_vel(n, vel_y_ant, vel_y_post, paso, fuerza_y_ant, fuerza_y_post);
    adv_vel(n, vel_z_ant, vel_z_post, paso, fuerza_z_ant, fuerza_z_post);
  }
  return 0;
}


// calcula la temperatura a partir del valor medio de la energia cinetica
// y tambien la energia cinetica, y guarda ambas en pointers
int kinetic_energy(int n_part, int iter, float* vel_x, float* vel_y, float* vel_z,
                        float* kinetic){

    float KE = 0;
    // calcula la energia cinetica total (KE)
    for (int i = 0; i < n_part; i++) {
        KE += (vel_x[i]*vel_x[i] + vel_y[i]*vel_y[i] + vel_z[i]*vel_z[i]);
    }
    KE = KE/2;

    // temp[iter] = (2./3.)*KE/n_part;
    kinetic[iter] = KE;

    return 0;
}

float temperature(int n_part, float* vel_x, float* vel_y, float* vel_z){
    float T;
    float KE = 0;
    for (int i = 0; i < n_part; i++) {
        KE += (vel_x[i]*vel_x[i] + vel_y[i]*vel_y[i] + vel_z[i]*vel_z[i]);
    }
    KE /= 2;
    T = (2./3.)*KE/n_part;
    return T;
}

// calcula la energía potencial para todas las partículas
// le agrego el shift para suavizar el cut_off, ver Haile ss 5.1.2
int potential_energy(int n_part, int iter, float l, float* x, float* y, float* z, float r_cut2, float* potential){

  int i, j;
  float dx, dy, dz;
  float invr2, distancia2;
  float U = 0;
  for (i = 0; i < n_part-1; i++) {
    for (j = i+1; j < n_part; j++) {
      dx = delta_coord(i, j, x, l);
      dy = delta_coord(i, j, y, l);
      dz = delta_coord(i, j, z, l);
      distancia2 = dist2(dx, dy, dz);
      invr2 = 1./distancia2;
      if (distancia2 < r_cut2) {
        U += 4 * pow(invr2,3) * ( pow(invr2,3) - 1. );
        U -= 4 * pow(1/r_cut2,3) * ( pow(1/r_cut2,3) - 1. );
		    //U += 4 * pow(invr2,3) * ( pow(invr2,3) - 1. ) + 0.0163 - (sqrt(distancia2)-2.5)*(-0.015599);
      }
    }
  }
  potential[iter] = U;

  return 0;
}

// calcula la energia total a partir de las otras
int total_energy(int iter, float* kinetic, float* potential, float* total){

  total[iter] = kinetic[iter] + potential[iter];

  return 0;
}

// hace el rescaling de velocidades para cambiar la temperatura
int rescaling(float T_inicial, float T_final, int n_part, float* vel_x, float* vel_y, float* vel_z){
	float scale = sqrt(T_final/T_inicial);
	int i;
	for (i = 0; i < n_part; i++){
		vel_x[i] *= scale;
		vel_y[i] *= scale;
		vel_z[i] *= scale;
	}

	// por las dudas les saca el valor medio..
	float s_x = 0;
	float s_y = 0;
	float s_z = 0;
	for (int n = 0; n < n_part; n++){
		s_x += vel_x[n];
		s_y += vel_y[n];
		s_z += vel_z[n];
	}
	for (int n = 0; n < n_part; n++){
		vel_x[n] -= s_x/n_part;
		vel_y[n] -= s_y/n_part;
		vel_z[n] -= s_z/n_part;
	}

	return 0;
}


int delta_N(float l, float r, float delta_r, int n_part, float* x, float* y, float* z){
  int i;
  float distancia2;
  int n = 0;
  for (i = 0; i < n_part; i++) {
    distancia2 = dist2(x[i] - l/2, y[i] - l/2, z[i] - l/2);
    if (distancia2 < pow(r + delta_r, 2) && distancia2 > pow(r, 2)) {
        n++;
    }
  }
  return n;
}

float delta_vol(float r, float delta_r){
  float volumen;
  volumen = (4./3.) * 3.14 * (pow(r + delta_r, 3) - pow(r, 3));
  return volumen;
}

float g(float rho, float l, float r, float delta_r, int n_part, float* x, float* y, float* z){
  float volumen = delta_vol(r, delta_r);
  int N = delta_N( l, r, delta_r, n_part, x, y, z);
  float g = N/(0.5 * n_part * rho * volumen);
  // printf("g = %f\n", g);
  return g;
}

float g_posta(float delta_r, int particiones_r, float l, float rho, int n_part, float* r, float* g, float* x, float* y, float* z){
	int i, j;
	float dx, dy, dz, distancia;
	int k;
	/*
	Recorre todos los pares de partículas y calcula su distancia. Expresa ésta como
	un múltiplo del dr, y le suma dos a la g(r) que corresponde a dicha distancia
	(por la contribución del par de partículas)
	*/
	for(i = 0; i < n_part-1; i++){
		for(j = i+1; j < n_part; j++){
			dx = delta_coord(i, j, x, l);
			dy = delta_coord(i, j, y, l);
			dz = delta_coord(i, j, z, l);
			distancia = sqrt(dist2(dx, dy, dz));
			k = (int) (distancia/delta_r);
			if (distancia < 0.5*l){
			g[k] += 2;
			}

		}
	}
	// ahora para obtener la g de verdad, divido esa cant de particulas por el volumen de la cascara y etc...
	float volumen;
	for(int k = 0; k < particiones_r; k++){
		volumen = delta_vol(r[k], delta_r);
		g[k] /= (0.5 * n_part * rho * volumen);
	}

}

float pressure(int n_part, float l, float rho, float* x, float* y, float* z, float r_cut2, float temperature){
  int i, j;
  float dx, dy, dz;
  float invr2, distancia2;
  float w = 0;
  float P = 0;
  for (i = 0; i < n_part-1; i++) {
    for (j = i+1; j < n_part; j++) {
      dx = delta_coord(i, j, x, l);
      dy = delta_coord(i, j, y, l);
      dz = delta_coord(i, j, z, l);
      distancia2 = dist2(dx, dy, dz);
      invr2 = 1./distancia2;
      if (distancia2 < r_cut2) {
          w += 48*(pow(invr2,3) - 0.5)*pow(invr2,3);
      }
    }
  }
  P = rho*(temperature + w/(3*n_part));
  return P;
}
