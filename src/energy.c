#include<stdlib.h>
#include<math.h>

#include "vector.h"
#include "energy.h"

#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))

void computeparticleEnergy2Particles(int particle1, int particle2, vector *pos, int numparticles, vector box, double *energy, double *virial){
	//private, only used in computeparticleEnergy and computeSystemEnergy
	double cutoff = CUTOFF;

	double dist2;
	vector dist;

	dist.x = pos[particle1].x - pos[particle2].x;
	dist.y = pos[particle1].y - pos[particle2].y;
	dist.z = pos[particle1].z - pos[particle2].z;

	dist.x -= box.x*rint(dist.x / box.x);
	dist.y -= box.y*rint(dist.y / box.y);
	dist.z -= box.z*rint(dist.z / box.z);

	dist2 = SQR(dist.x) + SQR(dist.y) + SQR(dist.z);

	if (dist2<SQR(cutoff))
	{
		double r2i = 1. / dist2;
		double r6i = CUBE(r2i);

		*energy += 4.*r6i*(r6i - 1.);
		*virial += 48.*r6i*(r6i - 0.5);
	}
}

int computeNumNeighbours(int particle1, int particle2, vector *pos, int numparticles, vector box, double *energy, double *virial){
	
	double cutoff = CUTOFF;

	double dist2;
	vector dist;

	dist.x = pos[particle1].x - pos[particle2].x;
	dist.y = pos[particle1].y - pos[particle2].y;
	dist.z = pos[particle1].z - pos[particle2].z;

	dist.x -= box.x*rint(dist.x / box.x);
	dist.y -= box.y*rint(dist.y / box.y);
	dist.z -= box.z*rint(dist.z / box.z);

	dist2 = SQR(dist.x) + SQR(dist.y) + SQR(dist.z);
	
	if (dist2<SQR(DIAMETER))
		return 1;
	else if (dist2<SQR(cutoff))
	{
		*energy = 1;
		//*virial += 48.*r6i*(r6i - 0.5);
	}
	else{
		*energy = 0;
	}
	return 0;
}

