#ifndef ENERGY_H
#define ENERGY_H

#define CUTOFF 1.5
#define DIAMETER 1.

void computeparticleEnergy2Particles(int particle1, int particle2, vector *pos, int numparticles, vector box, double *energy, double *virial);

int computeNumNeighbours(int particle1, int particle2, vector *pos, int numparticles, vector box, double *energy, double *virial);

#endif
