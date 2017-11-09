/*
 * sine_wave.c
 *
 *  Created on: Oct 23, 2017
 *      Author: docampo
 */





#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define PI  3.1415926535897931159979635
#define NP 500

int main (int argc, char *argv[])
{
	if (argc != 2) {
		printf(" Usage is %s + #periods (integers)\n",argv[0]);
		exit(1);
	}
	int n = atoi(argv[1]);
	
	int i;
	double tx, ty, tz;
	FILE *fil  = fopen("sineWave.dat","w");
	fprintf(fil,"#N: %d\n",NP);
	for ( i = 0 ; i < NP ; ++i) {
		tx = 2.0*(double)i/((double)NP-1.0);
		ty = sin((double)n *tx* PI);// + sin((double)m *ty* PI);
		tz = cos((double)n *tx* PI);// + sin((double)m *ty* PI);
		fprintf(fil,"%lf\t%lf\t%lf\n",tx,ty,tz);
      	}
	fclose(fil);
	return 0;
}

