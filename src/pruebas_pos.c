#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "pruebas_pos.h"




int init_prueba_dos(float l, float* pos_x_ant, float* pos_y_ant, float* pos_z_ant){
	float a = l/3.;
	pos_x_ant[0] = a;
	pos_y_ant[0] = a;
	pos_z_ant[0] = a;
	pos_x_ant[1] = 2*a;
	pos_y_ant[1] = 2*a;
	pos_z_ant[1] = 2*a;
}

int init_prueba_tres(float l, float* pos_x_ant, float* pos_y_ant, float* pos_z_ant){
	float a = l/4.;
	pos_x_ant[0] = a;
	pos_y_ant[0] = a;
	pos_z_ant[0] = a;
	pos_x_ant[1] = 2*a;
	pos_y_ant[1] = 2*a;
	pos_z_ant[1] = 2*a;
	pos_x_ant[2] = 3*a;
	pos_y_ant[2] = 3*a;
	pos_z_ant[2] = 3*a;
}

int init_prueba_cuatro(float l, float* pos_x_ant, float* pos_y_ant, float* pos_z_ant){
	float a = l/5.;
	pos_x_ant[0] = a;
	pos_y_ant[0] = a;
	pos_z_ant[0] = a;
	pos_x_ant[1] = 2*a;
	pos_y_ant[1] = 2*a;
	pos_z_ant[1] = 2*a;
	pos_x_ant[2] = 3*a;
	pos_y_ant[2] = 3*a;
	pos_z_ant[2] = 3*a;
	pos_x_ant[3] = 4*a;
	pos_y_ant[3] = 4*a;
	pos_z_ant[3] = 4*a;
}

int init_prueba_cinco(float l, float* pos_x_ant, float* pos_y_ant, float* pos_z_ant){
	float a = l/5.;
	pos_x_ant[0] = a;
	pos_y_ant[0] = a;
	pos_z_ant[0] = a;
	pos_x_ant[1] = 2*a;
	pos_y_ant[1] = 2*a;
	pos_z_ant[1] = 2*a;
	pos_x_ant[2] = 3*a;
	pos_y_ant[2] = 3*a;
	pos_z_ant[2] = 3*a;
	pos_x_ant[3] = 4*a;
	pos_y_ant[3] = 4*a;
	pos_z_ant[3] = 4*a;
	pos_x_ant[4] = a;
	pos_y_ant[4] = a;
	pos_z_ant[4] = 4*a;
}

int init_prueba_seis(float l, float* pos_x_ant, float* pos_y_ant, float* pos_z_ant){
	float a = l/5.;
	pos_x_ant[0] = a;
	pos_y_ant[0] = a;
	pos_z_ant[0] = a;
	pos_x_ant[1] = 2*a;
	pos_y_ant[1] = 2*a;
	pos_z_ant[1] = 2*a;
	pos_x_ant[2] = 3*a;
	pos_y_ant[2] = 3*a;
	pos_z_ant[2] = 3*a;
	pos_x_ant[3] = 4*a;
	pos_y_ant[3] = 4*a;
	pos_z_ant[3] = 4*a;
	pos_x_ant[4] = a;
	pos_y_ant[4] = a;
	pos_z_ant[4] = 4*a;
	pos_x_ant[5] = 2*a;
	pos_y_ant[5] = 2*a;
	pos_z_ant[5] = 4*a;
}