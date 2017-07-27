#include "stdlib.h"
#include <time.h>
#include <stdio.h>
#include "dynamics.h"
#include "pruebas_pos.h"
#include <math.h>

#define R_CUT2 6.25   // el posta es 6.25


int main(){

	FILE *fdat;
	fdat = fopen("prueba_fuerza_part.csv", "w");
	fprintf(fdat, "r	f(r)	U\n");
	
	float x[2];
	x[0] = 0.0;
	float F, U;
	float paso = 0.002;
	float l = 3;
	int iteraciones = 2.5/paso;
	float dx, distancia2, invr2;
	
	for (int i = 1; i < iteraciones; i++){
		
		x[1] = paso*i;
		dx = delta_coord(0, 1, x, l);
		distancia2 = dist2(dx, 0, 0);
		F = eval_LJ(distancia2, R_CUT2)*dx;
		invr2 = 1/distancia2;
		U = 4 * pow(invr2,3) * ( pow(invr2,3) - 1. );
		fprintf(fdat,"%.4g,%.4g,%.4g\n", dx, F, U);
		
	}
	
 
 
	return 0;
}