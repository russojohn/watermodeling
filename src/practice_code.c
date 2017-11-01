// add ;_CRT_SECURE_NO_WARNINGS to pre proccessor to dissable warnings

//#include "stdafx.h"
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<assert.h>
#include<string.h>

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
	int overlap=0;
	
	int particle1;
	
	*energy=0;
	
	for (particle1 = 0; particle1<numparticles; particle1++)
	{
		int neighbours=0;
		
		int particle2;
		for (particle2 = 0; particle2<numparticles; particle2++)
		{
			if (particle1!=particle2)
			{
			
				overlap+=computeNumNeighbours(particle1, particle2, pos, numparticles, box, &e, &v);
				
				neighbours+=(int)e;
			}
		}
		
		if (neighbours==4)
		{
			*energy+=-1;
		}
	}
	
	
	if (overlap>0)
	{
		fprintf(stderr,"Error in initial condition. There are %d overlap\n",overlap);
		exit(1);
	}
	
	
	
// 	*energy = e;
// 	*virial = v;

}

void computeparticleEnergy(int seed, vector *pos, int numparticles, vector box, double *energy, double *virial)
{
	//e and v are needed to prevent overflow/underflow
	double e = 0, v = 0;

	int particle;
	for (particle = 0; particle < seed; particle++)
	{
		computeparticleEnergy2Particles(seed, particle, pos, numparticles, box, &e, &v);
	}

	for (particle = seed + 1; particle < numparticles; particle++)
	{
		computeparticleEnergy2Particles(seed, particle, pos, numparticles, box, &e, &v);
	}

	*energy = e;
	*virial = v;
}


int computeparticleEnergyCells(int seed, vector *pos, int numparticles, vector box, double *energy, double *virial,listcell *cells)
{
	double e = 0, v = 0;
	int overlap=0;
	
	overlap=getParticleEnergy(cells,pos,numparticles,seed,box,&e);
	
	
	
	*energy = e;
	*virial = v;
	
	return overlap;
}


void metropolisMove(int seed, vector *pos, int numparticles, vector box, double maxdisplacement, double temperature, double *energy, double *virial,listcell *cells)
{
	double old_energy, old_virial;
	int overlap=0;
 	//computeparticleEnergyCells(seed, pos, numparticles, box, &old_energy, &old_virial,cells);
	overlap=getParticleEnergy(cells,pos,numparticles,seed,box,&old_energy);
	//computeparticleEnergy(seed, pos, numparticles, box, &old_energy, &old_virial);

	vector old_pos = pos[seed];

	pos[seed].x += maxdisplacement*(((double)rand() / (double)RAND_MAX) - 0.5);
	pos[seed].y += maxdisplacement*(((double)rand() / (double)RAND_MAX) - 0.5);
	pos[seed].z += maxdisplacement*(((double)rand() / (double)RAND_MAX) - 0.5);

	
	
	//updateList(cells,pos,numparticles);
	
	double new_energy, new_virial;

 	//computeparticleEnergyCells(seed, pos, numparticles, box, &new_energy, &new_virial,cells);
	
	overlap=getParticleEnergy(cells,pos,numparticles,seed,box,&new_energy);
	
	
	
	//computeparticleEnergy(seed, pos, numparticles, box, &new_energy, &new_virial);

	// metropolis
	double deltaE = new_energy - old_energy;
	
	if (overlap!=0)
		pos[seed] = old_pos;
	else
	{

		if ((deltaE <= 0.) || (((double)rand() / (double)RAND_MAX)<exp(-(deltaE) / temperature)))
		{
	// 		printf("%lf %lf %lf %lf\n",old_energy,new_energy,deltaE,exp(-(deltaE) / temperature));
			
			// accepted move
			*energy += deltaE;
			*virial += new_virial - old_virial;
			
			changeCell(cells,&old_pos,pos+seed,seed);
		}
		else
		{
			
			
			pos[seed] = old_pos;
		}
	}
}

void printStatistics(int t, vector *pos, int numparticles, vector box, double temperature, double energy, double virial, FILE *efile, FILE *pfile)
{
	double density = numparticles / (box.x*box.y*box.z);

	double ideal_pressure = density*temperature;

	double pressure = ideal_pressure + (virial) / (3.*box.x*box.y*box.z);

	double system_energy;
	computeSystemEnergy(pos, numparticles, box, &system_energy, &virial);
	
	fprintf(efile, "%d %9.5lf\n", t, (system_energy) / (double)numparticles);
	fprintf(pfile, "%d %9.5lf\n", t, pressure);
}

int main(int argc, char *argv[])
{
	if (argc!=2)
	{
		printf("%s [initial file]\n",argv[0]);
		exit(1);
	}
	
	// DEFINE VARIABLES
	char input_file[200];
	
	strcpy(input_file,argv[1]);
	
	vector *pos;
	int numparticles;
	vector box;
	int i; //why am i defined up here?
	int total_length = 10000000;
	double maxdisplacement = 0.1;
	double temperature = 0.5;

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
		int s;
		for (s = 0; s < numparticles; s++)
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
			
			
// 			
			//why dont we flush pfile?
		}
		
// 		printf("pos %lf %lf %lf\n",pos[99].x,pos[99].y,pos[99].z);
		
	}

	fclose(efile);
	fclose(pfile);

	return 0;
}