// add ;_CRT_SECURE_NO_WARNINGS to pre proccessor to dissable warnings

#include "stdafx.h"
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<assert.h>

#include "vector.h"
#include "practice_energy.h"
#include "practice_cells.h"

#define MAX_LINE_LENGTH 1000
#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))

void readInitialConditions(char *input_file, vector **pos, int *numparticles, vector *box, int *time)
{
	FILE *ifile = fopen(input_file, "r");

	if (ifile == NULL)
	{
		//why use fprintf?
		fprintf(stderr, "Error: can't open initial conditions file '%s'\n", input_file);
		exit(1);
	}


	// READ HEADER
	//NOTES: lf = long float

	char line[MAX_LINE_LENGTH] = "";
	fgets(line, MAX_LINE_LENGTH, ifile);

	int state = sscanf(line, "%d %d %lf %lf %lf\n", time, numparticles, &box->x, &box->y, &box->z);

	if (state != 5)
	{
		fprintf(stderr, "Error while reading header in '%s' file\n", input_file);
		exit(1);
	}


	// ALLOCATE MEMORY
	//*pos = new vector[*numparticles];
	*pos=calloc(*numparticles,sizeof(vector));


	// READ POSITIONS
	int n = 0;
	while (fgets(line, MAX_LINE_LENGTH, ifile))
	{
		state = sscanf(line, "%lf %lf %lf\n", &((*pos)[n].x), &((*pos)[n].y), &((*pos)[n].z));
		n++;

		if (state != 3)
		{
			fprintf(stderr, "Error while reading '%s' file\n", input_file);
			exit(1);
		}
	}

	assert(n == *numparticles);

	fclose(ifile);
}


void computeSystemEnergy(vector *pos, int numparticles, vector box, double *energy, double *virial)
{
	//no need to optimise, only called once, be lazy :)
	double e = 0, v = 0;

	for (int particle1 = 0; particle1<numparticles; particle1++)
	{
		for (int particle2 = particle1 + 1; particle2<numparticles; particle2++)
		{
			computeparticleEnergy2Particles(particle1, particle2, pos, numparticles, box, &e, &v);
		}
	}

	*energy = e;
	*virial = v;

}

void computeparticleEnergy(int seed, vector *pos, int numparticles, vector box, double *energy, double *virial)
{
	//e and v are needed to prevent overflow/underflow
	double e = 0, v = 0;

	for (int particle = 0; particle < seed; particle++)
	{
		computeparticleEnergy2Particles(seed, particle, pos, numparticles, box, &e, &v);
	}

	for (int particle = seed + 1; particle < numparticles; particle++)
	{
		computeparticleEnergy2Particles(seed, particle, pos, numparticles, box, &e, &v);
	}

	*energy = e;
	*virial = v;
}


void computeparticleEnergyCells(int seed, vector *pos, int numparticles, vector box, double *energy, double *virial,listcell *cells)
{
	double e = 0, v = 0;
	
	getParticleEnergy(cells,pos,numparticles,seed,box,&e);
	
}

void metropolisMove(int seed, vector *pos, int numparticles, vector box, double maxdisplacement, double temperature, double *energy, double *virial,listcell *cells)
{
	double old_energy, old_virial;

	computeparticleEnergyCells(seed, pos, numparticles, box, &old_energy, &old_virial,cells);

	vector old_pos = pos[seed];

	pos[seed].x += maxdisplacement*(((double)rand() / (double)RAND_MAX) - 0.5);
	pos[seed].y += maxdisplacement*(((double)rand() / (double)RAND_MAX) - 0.5);
	pos[seed].z += maxdisplacement*(((double)rand() / (double)RAND_MAX) - 0.5);

	double new_energy, new_virial;

	computeparticleEnergyCells(seed, pos, numparticles, box, &new_energy, &new_virial,cells);

	// metropolis
	double deltaE = new_energy - old_energy;

	if ((deltaE <= 0.) || (((double)rand() / (double)RAND_MAX)<exp(-(deltaE) / temperature)))
	{
		// accepted move
		*energy += deltaE;
		*virial += new_virial - old_virial;
	}
	else
	{
		pos[seed] = old_pos;
	}
}

void printStatistics(int t, vector *pos, int numparticles, vector box, double temperature, double energy, double virial, FILE *efile, FILE *pfile)
{
	double density = numparticles / (box.x*box.y*box.z);

	double ideal_pressure = density*temperature;

	double pressure = ideal_pressure + (virial) / (3.*box.x*box.y*box.z);

	fprintf(efile, "%d %9.5lf\n", t, (energy) / (double)numparticles);
	fprintf(pfile, "%d %9.5lf\n", t, pressure);
}

int main(int argc, char *argv[])
{
	// DEFINE VARIABLES
	char input_file[200] = "initial_conditions.dat";
	vector *pos;
	int numparticles;
	vector box;
	int i; //why am i defined up here?
	int total_length = 100000;
	double maxdisplacement = 0.1;
	double temperature = 1.2;

	// READ INITIAL FILE
	//why do we pass i? what is i?
	readInitialConditions(input_file, &pos, &numparticles, &box, &i);
	
	
	listcell *cells=getList(box,CUTOFF,numparticles);
	updateList(cells,pos,numparticles);
	
	
	// COMPUTE INITIAL ENERGY
	double energy, virial;
	computeSystemEnergy(pos, numparticles, box, &energy, &virial);
	
	// OPEN ENERGY AND PRESSURE FILES
	FILE *efile = fopen("energy.dat", "w");
	FILE *pfile = fopen("pressure.dat", "w");

	// INITIALIZE RANDOM NUMBERS
	time_t t;
	srand((unsigned)time(&t));

	while (++i<total_length)
	{
		// MONTE CARLO SWEEP
		for (int s = 0; s < numparticles; s++)
		{
			int seed = rand() % numparticles;
			metropolisMove(seed, pos, numparticles, box, maxdisplacement, temperature, &energy, &virial,cells);
		}

		if (i % 100 == 0)
		{
			printf("Step %d\n", i);
			printStatistics(i, pos, numparticles, box, temperature, energy, virial, efile, pfile);
			fflush(efile);
			fflush(pfile);
			//why dont we flush pfile?
		}
	}

	fclose(efile);
	fclose(pfile);

	return 0;
}