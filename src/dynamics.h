#ifndef DYNAMICS_H
#define DYNAMICS_H
int hola();
int initizalize(int n_part, float* pos_x_ant, float* pos_y_ant, float* pos_z_ant,
                float* vel_x_ant, float* vel_y_ant, float* vel_z_ant);
int adv_pos(float* pos_ant, float* pos_post, float* vel, float paso, float* fuerza);
int adv_vel(float* vel_ant, float* vel_post, float paso, float* fuerza_ant, float* fuerza_post);
#endif
