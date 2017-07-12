#ifndef DYNAMICS_H
#define DYNAMICS_H
int hola();
int init_prueba_dos(float l, float* pos_x_ant, float* pos_y_ant, float* pos_z_ant);
int initizalize_pos(int n_part, float l, float* pos_x_ant, float* pos_y_ant, float* pos_z_ant);
float gauss(float T);
int initizalize_vel(int n_part, float* vel_x_ant, float* vel_y_ant, float* vel_z_ant, float T);
int sum_vel(int n_part, float* vel_x, float* vel_y, float* vel_z);
int adv_pos(int n, float* pos_ant, float* pos_post, float l, float* vel, float paso, float paso2, float* fuerza);
int adv_vel(int n, float* vel_ant, float* vel_post, float paso, float* fuerza_ant, float* fuerza_post);
int sign(float a, float b);
float delta_coord(int i, int j, float* coord, float l);
float dist2(float dx, float dy, float dz);
float eval_LJ(float dist2, float r_cut);
int F_todas(int n_part, float l, float* x, float* y, float* z, float r_cut2, float *fuerza_x, float *fuerza_y, float *fuerza_z);
int F_tot(int i, int n_part, float l, float* x, float* y, float* z, float r_cut2, float *fuerza_x, float *fuerza_y, float *fuerza_z);
int time_evol(int n_part, float l, float paso, float paso2, float* pos_x_ant, float* pos_x_post,
              float* pos_y_ant, float* pos_y_post, float* pos_z_ant,
              float* pos_z_post, float* vel_x_ant, float* vel_x_post,
              float* vel_y_ant, float* vel_y_post, float* vel_z_ant,
              float* vel_z_post, float* fuerza_x_ant, float* fuerza_x_post,
              float* fuerza_y_ant, float* fuerza_y_post, float* fuerza_z_ant,
              float* fuerza_z_post, float r_cut2);
int kinetic_temperature(int n_part, int iter, float* vel_x, float* vel_y, float* vel_z, float* temp, float* kinetic);
int potential_energy(int n_part, int iter, float l, float* x, float* y, float* z, float r_cut2, float* potential);
int total_energy(int iter, float* kinetic, float* potential, float* total);
#endif
