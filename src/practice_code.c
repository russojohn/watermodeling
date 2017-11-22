// add ;_CRT_SECURE_NO_WARNINGS to pre proccessor to dissable warnings

//#include "stdafx.h"
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<assert.h>
#include<string.h>

#include "vector.h"
#include "linked_list.h"
#include "practice_energy.h"
#include "practice_cells.h"

#define MAX_LINE_LENGTH 1000
#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))

static double Total_energy;

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
	*pos = calloc(*numparticles, sizeof(vector));


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
	int overlap = 0;

	int particle1;

	*energy = 0;

	for (particle1 = 0; particle1<numparticles; particle1++)
	{
		int neighbours = 0;

		int particle2;
		for (particle2 = 0; particle2<numparticles; particle2++)
		{
			if (particle1 != particle2)
			{
				overlap += computeNumNeighbours(particle1, particle2, pos, numparticles, box, &e, &v);

				neighbours += (int)e;
			}
		}

		if (neighbours == 4)
		{
			*energy += -1;
		}
	}

	if (overlap>0)
	{
		fprintf(stderr, "Error in initial condition. There are %d overlap\n", overlap);
		exit(1);
	}

	// 	*energy = e;
	// 	*virial = v;

}


void metropolisMove(int seed, vector *pos, double* numNeighbours, int numparticles, vector box, double maxdisplacement, double temperature, double *energy, double *virial, listcell *cells)
{
	// list1/2 keep track of the neighbours of seed before and after the move
	node *list1 = 0; //dec
	node *list2 = 0; //inc

	int overlap = 0;

	double old_energy, old_virial = 0;
	overlap = getParticleEnergy(cells, pos, numparticles, seed, box, &old_energy,&list1);

	vector old_pos = pos[seed];

	pos[seed].x += maxdisplacement*(((double)rand() / (double)RAND_MAX) - 0.5);
	pos[seed].y += maxdisplacement*(((double)rand() / (double)RAND_MAX) - 0.5);
	pos[seed].z += maxdisplacement*(((double)rand() / (double)RAND_MAX) - 0.5);

	double new_energy, new_virial = 0;
	overlap = getParticleEnergy(cells, pos, numparticles, seed, box, &new_energy,&list2);

	//print_list(list1);
	//print_list(list2);

	//remove duplicates
	node *curr1 = list1, *curr2 = 0;

	JUMP:;
	while (curr1 != 0){
		curr2 = list2;

		while (curr2 != 0){
			if ((curr1->data) == (curr2->data)){
				//remove duplicate and advance curr1, curr2 through linked list
				int temp = curr1->data;

				remove_list(&list1, temp);
				remove_list(&list2, temp);

				curr1 = list1;
				curr2 = list2;

				//print_list(list1);
				//print_list(list2);
				//break; //can i use lable to continue at the first loop like in java?
				//use of go to to mimic lable then break in java, fix later?
				goto JUMP;
			}
			else{
				curr2 = curr2->next;
			}
		}
		
		curr1 = curr1->next;
	}

	//print_list(list1);
	//print_list(list2);

	//calculate energy
	//note, old_energy is now number of nei of particle seed before moving and new_energy after move
	int oldE=0, newE = 0;

	if (old_energy == 4)
		oldE -= 1;
	if (new_energy == 4)
		newE -= 1;

	curr1 = list1; 
	curr2 = list2;
	while (curr1 != 0)
	{
		if (numNeighbours[curr1->data] == 5){
			//numNeighbours[curr1->data] - 1 == 4
			newE -= 1;
		}
		if (numNeighbours[curr1->data] == 4){
			//numNeighbours[curr1->data] - 1 == 4
			newE += 1;
		}
		curr1 = curr1->next;
	}
	while (curr2 != 0)
	{
		if (numNeighbours[curr2->data] == 3){
			//numNeighbours[curr2->data] + 1 == 4
			newE -= 1;
		}
		if (numNeighbours[curr2->data] == 4){
			//numNeighbours[curr2->data] + 1 == 4
			newE += 1;
		}
		curr2 = curr2->next;
	}
	
	// metropolis
	double deltaE = (double)(newE-oldE);

	if (overlap != 0){
		pos[seed] = old_pos;
	}
	else{
		if ((deltaE <= 0.) || (((double)rand() / (double)RAND_MAX)<exp(-(deltaE) / temperature)))
		{
			// accepted move
			*energy += deltaE;
			*virial += new_virial - old_virial;

			changeCell(cells, &old_pos, pos + seed, seed);
			
			//update numNeighbours array
			numNeighbours[seed] = new_energy;

			curr1 = list1;
			curr2 = list2;
			while (curr1 != 0)
			{
				numNeighbours[curr1->data] -= 1;
				curr1 = curr1->next;
			}
			while (curr2 != 0)
			{
				numNeighbours[curr2->data] += 1;
				curr2 = curr2->next;
			}


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

	fprintf(efile, "%d %9.5lf %9.5lf\n", t, (system_energy) / (double)numparticles,energy/(double)numparticles);
	fprintf(pfile, "%d %9.5lf\n", t, pressure);
}

void initNumNeighbours(vector *pos, int numparticles, vector box, double* numNeighbours, listcell *l){

	double e = 0;
	int overlap = 0;

	node *list = 0, *curr = 0;

	int i;
	for (i = 0; i < numparticles; i++){
		overlap += getParticleEnergy(l, pos, numparticles, i, box, &e, &list);
		//make another getParticleEnergy without list?
		while (list != 0){
			remove_list(&list, list->data);
		}

		numNeighbours[i] = e;
		//add error for overlap
		//printf("%d\n", (int)numNeighbours[i]);
	}

	
};

int main(int argc, char *argv[])
{
	/*if (argc!=2)
	{
	printf("%s [initial file]\n",argv[0]);
	exit(1);
	}*/

	// DEFINE VARIABLES
	char input_file[200] = "initial_conditions1.dat";

	//strcpy(input_file,argv[1]);
 

	vector *pos;
	double* numNeighbours;
	int numparticles;
	vector box;
	int i;
	int total_length = 10000000;
	double maxdisplacement = 0.1;
	double temperature = 0.5;

	// READ INITIAL FILE
	readInitialConditions(input_file, &pos, &numparticles, &box, &i);

	// INIT LISTCELLS
	listcell *cells = getList(box, CUTOFF, numparticles);
	updateList(cells, pos, numparticles);

	// INIT NUM NEI
	numNeighbours = malloc(numparticles*sizeof(double));
	initNumNeighbours(pos, numparticles, box, numNeighbours, cells);

	// COMPUTE INITIAL ENERGY
	double energy, virial;
	computeSystemEnergy(pos, numparticles, box, &energy, &virial);

	Total_energy=energy;

	// OPEN ENERGY AND PRESSURE FILES
	FILE *efile = fopen("energy.dat", "w");
	FILE *pfile = fopen("pressure.dat", "w");

	// INITIALIZE RANDOM NUMBERS
	time_t t;
	srand((unsigned)time(&t));

	int s, step = 100;
	while (++i<total_length)
	{

		// MONTE CARLO SWEEP
		for (s = 0; s < numparticles; s++)
		{
			int seed = rand() % numparticles;
			metropolisMove(seed, pos,numNeighbours, numparticles, box, maxdisplacement, temperature, &energy, &virial, cells);
		}

		if (i % step == 0)
		{
			printf("Step %d\n", i);
			printStatistics(i, pos, numparticles, box, temperature, energy, virial, efile, pfile);
			fflush(efile);
			fflush(pfile);

		}
	}

	fclose(efile);
	fclose(pfile);

	return 0;
}
