/*
 * sine_wave.c
 *
 *  Created on: Oct 23, 2017
 *      Author: docampo
 */





#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define PI  3.1415926535897931159979635
//#define NP 500
#define HELIX      0
#define SPHERE_ARC 1
#define SPHERE_SNAKE 2
int main (int argc, char *argv[])
{
	if (argc != 3) {
		printf(" Usage is %s + SHAPE: 0 = HELIX or 1 = SPHERE_ARC 2= SPHERE_SNAKE + NUMBER OF POINTS \n",argv[0]);
		exit(1);
	}
	int shape = atoi(argv[1]);
	int NP    = atoi(argv[2]);
	int i, n;
	double tx, ty, tz, theta, phi;
	FILE *fil;
	switch(shape) {
	case HELIX:
	{
		fil = fopen("helix.dat","w");
		fprintf(fil,"#N: %d\n",NP);
		printf(" enter integer number: Periods for the helix\n");
		scanf("%d",&n);
		for ( i = 0 ; i < NP ; ++i) {
			tx = 2.0*(double)i/((double)NP-1.0);
			ty = sin((double)n *tx* PI);// + sin((double)m *ty* PI);
			tz = cos((double)n *tx* PI);// + sin((double)m *ty* PI);
			fprintf(fil,"%lf\t%lf\t%lf\n",tx,ty,tz);
		}
		fclose(fil);
		break;
	}
	case SPHERE_ARC:
	{
		fil = fopen("sphere_arc.dat","w");
		fprintf(fil,"#N: %d\n",NP);
		theta = M_PI*0.25;
		for ( i = 0 ; i < NP ; ++i) {
			phi = (double)i*M_PI/((double)NP-1.0);
			tx = cos(theta)*sin(phi);
			ty = sin(theta)*sin(phi);
			tz = cos(phi);
			fprintf(fil,"%lf\t%lf\t%lf\n",tx,ty,tz);
		}
		fclose(fil);
		break;
	}
	case SPHERE_SNAKE:
		{
			fil = fopen("sphere_snake.dat","w");
			fprintf(fil,"#N: %d\n",NP);
			theta = M_PI*0.25;
			for ( i = 0 ; i < NP ; ++i) {
				theta = (double)i*M_PI/((double)NP/4-1.0);
				phi = (double)i*M_PI/((double)NP-1.0);
				tx = cos(theta)*sin(phi);
				ty = sin(theta)*sin(phi);
				tz = cos(phi);
				fprintf(fil,"%lf\t%lf\t%lf\n",tx,ty,tz);
			}
			fclose(fil);
			break;
		}

	}
	return 0;
}

